#!/usr/bin/env python
#
# Copyright (c) 2020-2030 Translational Oncology at the Medical Center of the Johannes Gutenberg-University Mainz gGmbH.
#
# This file is part of Neofox
# (see https://github.com/tron-bioinformatics/neofox).
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program. If not, see <http://www.gnu.org/licenses/>.#
import logging
import os
import time
from typing import List
import logzero
from logzero import logger
from dask.distributed import Client

from neofox.published_features.Tcell_predictor.tcellpredictor_wrapper import TcellPrediction
from neofox.references.references import ReferenceFolder, DependenciesConfiguration
from neofox import NEOFOX_LOG_FILE_ENV
from neofox.annotator import NeoantigenAnnotator
from neofox.exceptions import NeofoxConfigurationException, NeofoxDataValidationException
from neofox.model.neoantigen import NeoantigenAnnotations, Neoantigen, Patient
from neofox.model.conversion import ModelValidator


class NeoFox:

    def __init__(self, neoantigens: List[Neoantigen], patient_id: str, patients: List[Patient], num_cpus: int,
                 work_folder=None, output_prefix=None, reference_folder: ReferenceFolder = None,
                 configuration: DependenciesConfiguration = None):

        # initialise logs
        self._initialise_logs(output_prefix, work_folder)

        # intialize references folder and configuration
        # NOTE: uses the reference folder and config passed as a parameter if exists, this is here to make it
        # testable with fake objects
        self.reference_folder = reference_folder if reference_folder else ReferenceFolder()
        self.configuration = configuration if configuration else DependenciesConfiguration()
        self.tcell_predictor = TcellPrediction()
        self.num_cpus = num_cpus

        if neoantigens is None or len(neoantigens) == 0 or patients is None or len(patients) == 0:
            raise NeofoxConfigurationException("Missing input data to run Neofox")

        # TODO: avoid overriding patient id parameter
        for n in neoantigens:
            if n.patient_identifier is None:
                n.patient_identifier = patient_id

        # validates input data
        self.neoantigens = [ModelValidator.validate_neoantigen(n) for n in neoantigens]
        self.patients = {patient.identifier: ModelValidator.validate_patient(patient) for patient in patients}
        self._validate_input_data()

        logger.info("Data loaded")

    def _initialise_logs(self, output_prefix, work_folder):
        if work_folder and os.path.exists(work_folder):
            logfile = os.path.join(work_folder, "{}.log".format(output_prefix))
        else:
            logfile = os.environ.get(NEOFOX_LOG_FILE_ENV)
        if logfile is not None:
            logzero.logfile(logfile)
        # TODO: this does not work
        logzero.loglevel(logging.DEBUG)
        logger.info("Loading data...")

    def _validate_input_data(self):

        patient_identifiers_from_neoantigens = set([n.patient_identifier for n in self.neoantigens])
        patient_identifiers_from_patients = set([p.identifier for p in self.patients.values()])

        # check that there are no repeated neoantigens
        neoantigen_identifiers = [n.identifier for n in self.neoantigens]
        if len(neoantigen_identifiers) != len(set(neoantigen_identifiers)):
            raise NeofoxDataValidationException("There are repeated neoantigens!")

        # checks that no neoantigen is referring to an empty patient
        if "" in patient_identifiers_from_neoantigens or None in patient_identifiers_from_neoantigens:
            raise NeofoxDataValidationException(
                "There are neoantigens missing a reference to a patient")

        # checks that there is no neoantigen referring to a non existing patient
        missing_patient_identifiers = patient_identifiers_from_neoantigens.difference(patient_identifiers_from_patients)
        if len(missing_patient_identifiers) > 0:
            raise NeofoxDataValidationException(
                "There are neoantigens referring to missing patients: {}".format(missing_patient_identifiers))

    def get_annotations(self) -> List[NeoantigenAnnotations]:
        """
        Loads epitope data (if file has been not imported to R; colnames need to be changed), adds data to class that are needed to calculate,
        calls epitope class --> determination of epitope properties,
        write to txt file
        """
        logger.info("Starting NeoFox annotations...")
        # initialise dask
        # TODO: number of threads is hard coded. Is there a better value for this?
        dask_client = Client(processes=True, n_workers=self.num_cpus, threads_per_worker=4)
        # feature calculation for each epitope
        futures = []
        start = time.time()
        for neoantigen in self.neoantigens:
            patient = self.patients.get(neoantigen.patient_identifier)
            logger.debug("Neoantigen: {}".format(neoantigen.to_json(indent=3)))
            logger.debug("Patient: {}".format(patient.to_json(indent=3)))
            futures.append(dask_client.submit(
                NeoFox.annotate_neoantigen, neoantigen, patient, self.reference_folder, self.configuration,
                self.tcell_predictor
            ))

        annotations = dask_client.gather(futures)
        dask_client.close()
        end = time.time()
        logger.info("Elapsed time for annotating {} neoantigens {} seconds".format(
            len(self.neoantigens), int(end - start)))
        return annotations

    @staticmethod
    def annotate_neoantigen(neoantigen: Neoantigen, patient: Patient, reference_folder: ReferenceFolder,
                            configuration: DependenciesConfiguration, tcell_predictor: TcellPrediction):
        logger.info("Starting neoantigen annotation: {}".format(neoantigen.identifier))
        start = time.time()
        annotation = NeoantigenAnnotator(reference_folder, configuration, tcell_predictor).get_annotation(neoantigen, patient)
        end = time.time()
        logger.info("Elapsed time for annotating neoantigen {}: {} seconds".format(
            neoantigen.identifier, int(end - start)))
        return annotation
