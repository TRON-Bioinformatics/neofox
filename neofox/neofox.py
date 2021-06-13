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
from dask.distributed import performance_report

from neofox.expression_imputation.expression_imputation import ExpressionAnnotator
from neofox.published_features.Tcell_predictor.tcellpredictor_wrapper import (
    TcellPrediction,
)
from neofox.published_features.self_similarity.self_similarity import (
    SelfSimilarityCalculator,
)
from neofox.references.references import ReferenceFolder, DependenciesConfiguration
from neofox import NEOFOX_LOG_FILE_ENV, AFFINITY_THRESHOLD_DEFAULT
from neofox.annotator import NeoantigenAnnotator
from neofox.exceptions import (
    NeofoxConfigurationException,
    NeofoxDataValidationException,
)
from neofox.model.neoantigen import NeoantigenAnnotations, Neoantigen, Patient
from neofox.model.conversion import ModelValidator
import dotenv


class NeoFox:

    def __init__(
        self,
        neoantigens: List[Neoantigen],
        patients: List[Patient],
        num_cpus: int = 1,
        patient_id: str = None,
        work_folder=None,
        output_prefix=None,
        reference_folder: ReferenceFolder = None,
        configuration: DependenciesConfiguration = None,
        verbose=False,
        configuration_file=None,
        affinity_threshold=AFFINITY_THRESHOLD_DEFAULT
    ):

        self.affinity_threshold = affinity_threshold

        if configuration_file:
            dotenv.load_dotenv(configuration_file, override=True)

        # initialise logs
        self.log_file_name = self._get_log_file_name(output_prefix, work_folder)
        self._initialise_logs(self.log_file_name, verbose)

        # intialize references folder and configuration
        # NOTE: uses the reference folder and config passed as a parameter if exists, this is here to make it
        # testable with fake objects
        self.reference_folder = (
            reference_folder if reference_folder else ReferenceFolder()
        )
        # NOTE: makes this call to force the loading of the available alleles here
        self.reference_folder.get_available_alleles()
        self.configuration = (
            configuration if configuration else DependenciesConfiguration()
        )
        self.tcell_predictor = TcellPrediction(affinity_threshold=self.affinity_threshold)
        self.self_similarity = SelfSimilarityCalculator()
        self.num_cpus = num_cpus

        if (
            neoantigens is None
            or len(neoantigens) == 0
            or patients is None
            or len(patients) == 0
        ):
            raise NeofoxConfigurationException("Missing input data to run Neofox")

        # TODO: avoid overriding patient id parameter
        for n in neoantigens:
            if n.patient_identifier is None:
                n.patient_identifier = patient_id

        # validates input data
        self.neoantigens = [ModelValidator.validate_neoantigen(n) for n in neoantigens]
        self.patients = {
            patient.identifier: ModelValidator.validate_patient(patient)
            for patient in patients
        }
        self._validate_input_data()

        # retrieve from the data, if RNA-seq was available
        # add this information to patient model
        expression_per_patient = {self.patients[patient].identifier: [] for patient in self.patients}
        for neoantigen in self.neoantigens:
            expression_per_patient[neoantigen.patient_identifier].append(neoantigen.rna_expression)

        for patient in self.patients:
            self.patients[patient].is_rna_available = all(e is not None for e in expression_per_patient[self.patients[patient].identifier])

        # impute expresssion from TCGA, ONLY if isRNAavailable = False for given patient,
        # otherwise original values is reported
        # NOTE: this must happen after validation to avoid uncaptured errors due to missing patients
        # NOTE: add gene expression to neoantigen candidate model
        self.neoantigens = self._conditional_expression_imputation()

        logger.info("Data loaded")

    def _conditional_expression_imputation(self) -> List[Neoantigen]:
        expression_annotator = ExpressionAnnotator()
        neoantigens_transformed = []

        for neoantigen in self.neoantigens:
            expression_value = neoantigen.rna_expression
            patient = self.patients[neoantigen.patient_identifier]
            neoantigen_transformed = neoantigen
            gene_expression = expression_annotator.get_gene_expression_annotation(
                gene_name=neoantigen.gene, tcga_cohort=patient.tumor_type
            )
            if not patient.is_rna_available and patient.tumor_type is not None and patient.tumor_type != "":
                expression_value = gene_expression
            neoantigen_transformed.rna_expression = expression_value
            neoantigen.imputed_gene_expression = gene_expression
            neoantigens_transformed.append(neoantigen_transformed)
        return neoantigens_transformed

    @staticmethod
    def _initialise_logs(logfile, verbose=False):
        if logfile is not None:
            logzero.logfile(logfile)
        # TODO: this does not work
        if verbose:
            logzero.loglevel(logging.DEBUG)
        else:
            logzero.loglevel(logging.INFO)

    def _get_log_file_name(self, output_prefix, work_folder):
        if work_folder and os.path.exists(work_folder):
            logfile = os.path.join(work_folder, "{}.log".format(output_prefix))
        else:
            logfile = os.environ.get(NEOFOX_LOG_FILE_ENV)
        return logfile


    def _validate_input_data(self):

        patient_identifiers_from_neoantigens = set(
            [n.patient_identifier for n in self.neoantigens]
        )
        patient_identifiers_from_patients = set(
            [p.identifier for p in self.patients.values()]
        )

        # check that there are no repeated neoantigens
        neoantigen_identifiers = [n.identifier for n in self.neoantigens]
        if len(neoantigen_identifiers) != len(set(neoantigen_identifiers)):
            raise NeofoxDataValidationException("There are repeated neoantigens!")

        # checks that no neoantigen is referring to an empty patient
        if (
            "" in patient_identifiers_from_neoantigens
            or None in patient_identifiers_from_neoantigens
        ):
            raise NeofoxDataValidationException(
                "There are neoantigens missing a reference to a patient"
            )

        # checks that there is no neoantigen referring to a non existing patient
        missing_patient_identifiers = patient_identifiers_from_neoantigens.difference(
            patient_identifiers_from_patients
        )
        if len(missing_patient_identifiers) > 0:
            raise NeofoxDataValidationException(
                "There are neoantigens referring to missing patients: {}".format(
                    missing_patient_identifiers
                )
            )

    def get_annotations(self) -> List[NeoantigenAnnotations]:
        """
        Loads epitope data (if file has been not imported to R; colnames need to be changed), adds data to class that are needed to calculate,
        calls epitope class --> determination of epitope properties,
        write to txt file
        """
        logger.info("Starting NeoFox annotations...")
        # initialise dask
        # see reference on using threads versus CPUs here https://docs.dask.org/en/latest/setup/single-machine.html
        dask_client = Client(
            n_workers=self.num_cpus, threads_per_worker=1,
        )
        annotations = self.send_to_client(dask_client)
        dask_client.close()

        return annotations

    def send_to_client(self, dask_client):
        # feature calculation for each epitope
        futures = []
        start = time.time()
        # NOTE: sets those heavy resources to be used by all workers in the cluster
        future_tcell_predictor = dask_client.scatter(
            self.tcell_predictor, broadcast=True
        )
        future_self_similarity = dask_client.scatter(
            self.self_similarity, broadcast=True
        )
        future_reference_folder = dask_client.scatter(
            self.reference_folder, broadcast=True
        )
        future_configuration = dask_client.scatter(self.configuration, broadcast=True)
        for neoantigen in self.neoantigens:
            patient = self.patients.get(neoantigen.patient_identifier)
            logger.debug("Neoantigen: {}".format(neoantigen.to_json(indent=3)))
            logger.debug("Patient: {}".format(patient.to_json(indent=3)))
            futures.append(
                dask_client.submit(
                    NeoFox.annotate_neoantigen,
                    neoantigen,
                    patient,
                    future_reference_folder,
                    future_configuration,
                    future_tcell_predictor,
                    future_self_similarity,
                    self.log_file_name,
                    self.affinity_threshold
                )
            )
        annotations = dask_client.gather(futures)
        end = time.time()
        logger.info(
            "Elapsed time for annotating {} neoantigens {} seconds".format(
                len(self.neoantigens), int(end - start)
            )
        )
        return annotations

    @staticmethod
    def annotate_neoantigen(
        neoantigen: Neoantigen,
        patient: Patient,
        reference_folder: ReferenceFolder,
        configuration: DependenciesConfiguration,
        tcell_predictor: TcellPrediction,
        self_similarity: SelfSimilarityCalculator,
        log_file_name: str,
        affinity_threshold = AFFINITY_THRESHOLD_DEFAULT
    ):
        # the logs need to be initialised inside every dask job
        NeoFox._initialise_logs(log_file_name)
        logger.info("Starting neoantigen annotation id='{}' and peptide={}".format(
            neoantigen.identifier, neoantigen.mutation.mutated_xmer)
        )
        start = time.time()
        try:
            annotation = NeoantigenAnnotator(
                reference_folder,
                configuration,
                tcell_predictor=tcell_predictor,
                self_similarity=self_similarity,
                affinity_threshold=affinity_threshold
            ).get_annotation(neoantigen, patient)
        except Exception as e:
            logger.error("Error processing neoantigen {}".format(neoantigen.to_dict()))
            logger.error("Error processing patient {}".format(patient.to_dict()))
            raise e
        end = time.time()
        logger.info(
            "Elapsed time for annotating neoantigen id='{}' and peptide={}: {} seconds".format(
                neoantigen.identifier, neoantigen.mutation.mutated_xmer, int(end - start))
        )
        return annotation
