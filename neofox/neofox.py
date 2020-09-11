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
from neofox import NEOFOX_LOG_FILE_ENV
from neofox.annotator import NeoantigenAnnotator
from neofox.exceptions import NeofoxConfigurationException
from neofox.model.neoantigen import NeoantigenAnnotations, Neoantigen, Patient


class NeoFox:


    def __init__(self, neoantigens: List[Neoantigen], patient_id: str, patients: List[Patient], num_cpus: int, work_folder=None,
                 output_prefix = None):

        # initialise logs
        if work_folder and os.path.exists(work_folder):
            logfile = os.path.join(work_folder, "{}.log".format(output_prefix))
        else:
            logfile = os.environ.get(NEOFOX_LOG_FILE_ENV)
        if logfile is not None:
            logzero.logfile(logfile)
        # TODO: this does not work
        logzero.loglevel(logging.DEBUG)
        logger.info("Loading data...")

        # initialise dask
        # TODO: number of threads is hard coded. Is there a better value for this?
        self.dask_client = Client(processes=True, n_workers=num_cpus, threads_per_worker=4)

        if neoantigens is None or patients is None:
            raise NeofoxConfigurationException("Missing input data to run Neofox")
        self.neoantigens = neoantigens
        self.patients = {patient.identifier: patient for patient in patients}
        # TODO: avoid overriding patient id parameter
        for n in self.neoantigens:
            if n.patient_identifier is None:
                n.patient_identifier = patient_id

        logger.info("Data loaded")

    def get_annotations(self) -> List[NeoantigenAnnotations]:
        """
        Loads epitope data (if file has been not imported to R; colnames need to be changed), adds data to class that are needed to calculate,
        calls epitope class --> determination of epitope properties,
        write to txt file
        """
        logger.info("Starting NeoFox annotations...")
        # feature calculation for each epitope
        futures = []
        start = time.time()
        for neoantigen in self.neoantigens:
            patient = self.patients.get(neoantigen.patient_identifier)
            logger.debug("Neoantigen: {}".format(neoantigen.to_json(indent=3)))
            logger.debug("Patient: {}".format(patient.to_json(indent=3)))
            futures.append(self.dask_client.submit(NeoFox.annotate_neoantigen, neoantigen, patient))

        annotations = self.dask_client.gather(futures)
        end = time.time()
        logger.info("Elapsed time for annotating {} neoantigens {} seconds".format(
            len(self.neoantigens), int(end - start)))
        return annotations

    @staticmethod
    def annotate_neoantigen(neoantigen: Neoantigen, patient: Patient):
        logger.info("Starting neoantigen annotation: {}".format(neoantigen.identifier))
        start = time.time()
        annotation = NeoantigenAnnotator().get_annotation(neoantigen, patient)
        end = time.time()
        logger.info("Elapsed time for annotating neoantigen {}: {} seconds".format(
            neoantigen.identifier, int(end - start)))
        return annotation
