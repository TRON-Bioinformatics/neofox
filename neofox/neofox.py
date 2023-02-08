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

import neofox
from neofox.expression_imputation.expression_imputation import ExpressionAnnotator
from neofox.model.factories import NeoantigenFactory
from neofox.published_features.Tcell_predictor.tcellpredictor_wrapper import TcellPrediction
from neofox.published_features.self_similarity.self_similarity import SelfSimilarityCalculator
from neofox.references.references import ReferenceFolder, DependenciesConfiguration, ORGANISM_HOMO_SAPIENS
from neofox import NEOFOX_LOG_FILE_ENV
from neofox.annotator.neoantigen_annotator import NeoantigenAnnotator
from neofox.exceptions import NeofoxConfigurationException, NeofoxDataValidationException
from neofox.model.neoantigen import Neoantigen, Patient
from neofox.model.validation import ModelValidator
import dotenv


class NeoFox:

    def __init__(
            self,
            neoantigens: List[Neoantigen],
            patients: List[Patient],
            num_cpus: int = 1,
            log_file_name=None,
            reference_folder: ReferenceFolder = None,
            configuration: DependenciesConfiguration = None,
            verbose=True,
            configuration_file=None,
            rank_mhci_threshold=neofox.RANK_MHCI_THRESHOLD_DEFAULT,
            rank_mhcii_threshold=neofox.RANK_MHCII_THRESHOLD_DEFAULT,
            with_all_neoepitopes=False):

        initialise_logs(logfile=log_file_name, verbose=verbose)
        logger.info("Loading reference data...")

        self.rank_mhci_threshold = rank_mhci_threshold
        self.rank_mhcii_threshold = rank_mhcii_threshold

        if configuration_file:
            dotenv.load_dotenv(configuration_file, override=True)

        self.log_file_name = log_file_name

        # intialize references folder and configuration
        # NOTE: uses the reference folder and config passed as a parameter if exists, this is here to make it
        # testable with fake objects
        self.reference_folder = (
            reference_folder if reference_folder else ReferenceFolder(verbose=verbose)
        )
        # NOTE: makes this call to force the loading of the available alleles here
        self.reference_folder.get_available_alleles()
        self.configuration = (
            configuration if configuration else DependenciesConfiguration()
        )
        self.tcell_predictor = TcellPrediction()
        self.self_similarity = SelfSimilarityCalculator()
        self.num_cpus = num_cpus

        if (
            neoantigens is None
            or len(neoantigens) == 0
            or patients is None
            or len(patients) == 0
        ):
            raise NeofoxConfigurationException("Missing input data to run Neofox")

        # validates neoantigens
        self.neoantigens = neoantigens
        for n in self.neoantigens:
            # NOTE: the position of the mutations is not expected from the user and if provide the value is ignored
            n.position = NeoantigenFactory.mut_position_xmer_seq(neoantigen=n)
            ModelValidator.validate_neoantigen(n)

        # validates patients
        self.patients = {}
        for patient in patients:
            ModelValidator.validate_patient(patient, organism=self.reference_folder.organism)
            self.patients[patient.identifier] = patient

        self._validate_input_data()

        # retrieve from the data, if RNA-seq was available
        # add this information to patient model
        expression_per_patient = {self.patients[patient].identifier: [] for patient in self.patients}
        for neoantigen in self.neoantigens:
            expression_per_patient[neoantigen.patient_identifier].append(neoantigen.rna_expression)

        # only performs the expression imputation for humans
        if self.reference_folder.organism == ORGANISM_HOMO_SAPIENS:
            # impute expresssion from TCGA, ONLY if isRNAavailable = False for given patient,
            # otherwise original values is reported
            # NOTE: this must happen after validation to avoid uncaptured errors due to missing patients
            # NOTE: add gene expression to neoantigen candidate model
            self.neoantigens = self._conditional_expression_imputation()

        self.with_all_neoepitopes = with_all_neoepitopes

        logger.info("Reference data loaded")

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
            if expression_value is None and patient.tumor_type is not None and patient.tumor_type != "":
                expression_value = gene_expression
            neoantigen_transformed.rna_expression = expression_value
            neoantigen.imputed_gene_expression = gene_expression
            neoantigens_transformed.append(neoantigen_transformed)
        return neoantigens_transformed

    @staticmethod
    def get_log_file_name(output_prefix, work_folder):
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

    def get_annotations(self) -> List[Neoantigen]:
        """
        Loads epitope data (if file has been not imported to R; colnames need to be changed), adds data to class that are needed to calculate,
        calls epitope class --> determination of epitope properties,
        write to txt file
        """
        logger.info("Starting NeoFox annotations...")
        # initialise dask
        # see reference on using threads versus CPUs here https://docs.dask.org/en/latest/setup/single-machine.html
        dask_client = Client(n_workers=self.num_cpus, threads_per_worker=1)
        annotations = self.send_to_client(dask_client)
        dask_client.shutdown()          # terminates schedulers and workers
        dask_client.close(timeout=10)   # waits 10 seconds for the client to close before killing

        return annotations

    def send_to_client(self, dask_client):
        # feature calculation for each epitope
        futures = []
        start = time.time()
        # NOTE: sets those heavy resources to be used by all workers in the cluster
        future_tcell_predictor = dask_client.scatter(
            self.tcell_predictor, broadcast=True
        )
        future_self_similarity = dask_client.scatter(self.self_similarity, broadcast=True)
        future_reference_folder = dask_client.scatter(self.reference_folder, broadcast=True)
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
                    self.rank_mhci_threshold,
                    self.rank_mhcii_threshold,
                    self.with_all_neoepitopes
                )
            )
        annotated_neoantigens = dask_client.gather(futures)
        end = time.time()
        logger.info(
            "Elapsed time for annotating {} neoantigens {} seconds".format(
                len(self.neoantigens), int(end - start)
            )
        )
        return annotated_neoantigens

    @staticmethod
    def annotate_neoantigen(
        neoantigen: Neoantigen,
        patient: Patient,
        reference_folder: ReferenceFolder,
        configuration: DependenciesConfiguration,
        tcell_predictor: TcellPrediction,
        self_similarity: SelfSimilarityCalculator,
        log_file_name: str,
        rank_mhci_threshold = neofox.RANK_MHCI_THRESHOLD_DEFAULT,
        rank_mhcii_threshold=neofox.RANK_MHCII_THRESHOLD_DEFAULT,
        with_all_neoepitopes=False
    ):
        # the logs need to be initialised inside every dask job
        initialise_logs(log_file_name)
        logger.info("Starting neoantigen annotation with peptide={}".format(neoantigen.mutated_xmer))
        start = time.time()
        try:
            annotated_neoantigen = NeoantigenAnnotator(
                reference_folder,
                configuration,
                tcell_predictor=tcell_predictor,
                self_similarity=self_similarity,
                rank_mhci_threshold=rank_mhci_threshold,
                rank_mhcii_threshold=rank_mhcii_threshold
            ).get_annotated_neoantigen(neoantigen, patient, with_all_neoepitopes=with_all_neoepitopes)
        except Exception as e:
            logger.error("Error processing neoantigen {}".format(neoantigen.to_dict()))
            logger.error("Error processing patient {}".format(patient.to_dict()))
            raise e
        end = time.time()
        logger.info(
            "Elapsed time for annotating neoantigen for peptide={}: {} seconds".format(
                neoantigen.mutated_xmer, int(end - start))
        )
        return annotated_neoantigen


def initialise_logs(logfile, verbose=False):
    if logfile is not None:
        logzero.logfile(logfile)
    # TODO: this does not work
    if verbose:
        logzero.loglevel(logging.INFO)
    else:
        logzero.loglevel(logging.WARN)
