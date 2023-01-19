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
import time
from copy import copy
from typing import List
import logzero
from logzero import logger
from dask.distributed import Client

from neofox.annotator.neoepitope_annotator import NeoepitopeAnnotator
from neofox.expression_imputation.expression_imputation import ExpressionAnnotator
from neofox.published_features.Tcell_predictor.tcellpredictor_wrapper import TcellPrediction
from neofox.published_features.self_similarity.self_similarity import SelfSimilarityCalculator
from neofox.references.references import ReferenceFolder, DependenciesConfiguration, ORGANISM_HOMO_SAPIENS
from neofox.exceptions import NeofoxConfigurationException, NeofoxDataValidationException
from neofox.model.neoantigen import Patient, PredictedEpitope
from neofox.model.validation import ModelValidator
import dotenv


class NeoFoxEpitope:

    def __init__(
            self,
            neoepitopes: List[PredictedEpitope],
            patients: List[Patient] = [],
            num_cpus: int = 1,
            log_file_name=None,
            reference_folder: ReferenceFolder = None,
            configuration: DependenciesConfiguration = None,
            verbose=True,
            configuration_file=None):

        initialise_logs(logfile=log_file_name, verbose=verbose)
        logger.info("Loading reference data...")

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

        # validates patients
        self.patients = {}
        for patient in patients:
            ModelValidator.validate_patient(patient, organism=self.reference_folder.organism)
            self.patients[patient.identifier] = patient

        if neoepitopes is None or len(neoepitopes) == 0:
            raise NeofoxConfigurationException("Missing input data to run Neofox")

        # validates neoepitopes and combines neoepitopes according to patient alleles
        self.neoepitopes = []
        for n in neoepitopes:
            ModelValidator.validate_neoepitope(n, organism=self.reference_folder.organism)
            if ModelValidator.is_mhci_epitope(n) or ModelValidator.is_mhcii_epitope(n):
                self.neoepitopes.append(n)
            else:
                if n.patient_identifier in self.patients:
                    patient = self.patients.get(n.patient_identifier)
                    if ModelValidator.is_mhci_peptide_length_valid(len(n.mutated_peptide)):
                         for m in patient.mhc1:
                             for a in m.alleles:
                                mhci_neoepitope = copy(n)
                                mhci_neoepitope.allele_mhc_i = a
                                self.neoepitopes.append(mhci_neoepitope)
                    for m in patient.mhc2:
                        for i in m.isoforms:
                            mhcii_neoepitope = copy(n)
                            mhcii_neoepitope.isoform_mhc_i_i = i
                            self.neoepitopes.append(mhcii_neoepitope)
                else:
                    raise NeofoxDataValidationException(
                        'A neoepitope is linked to patient {} for which there is no data'.format(n.patient_identifier))

        # only performs the expression imputation for humans
        if self.reference_folder.organism == ORGANISM_HOMO_SAPIENS:
            # impute expresssion from TCGA, ONLY if isRNAavailable = False for given patient,
            # otherwise original values is reported
            # NOTE: this must happen after validation to avoid uncaptured errors due to missing patients
            # NOTE: add gene expression to neoantigen candidate model
            self.neoepitopes = self._conditional_expression_imputation()

        logger.info("Reference data loaded")

    def get_annotations(self) -> List[PredictedEpitope]:
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
        dask_client.close(timeout=10)   # waits 10 seconds for the client to close before killing

        return annotations

    def send_to_client(self, dask_client):
        # feature calculation for each epitope
        futures = []
        start = time.time()
        # NOTE: sets those heavy resources distributed to all workers in the cluster
        future_tcell_predictor = dask_client.scatter(
            self.tcell_predictor, broadcast=True
        )
        future_self_similarity = dask_client.scatter(self.self_similarity, broadcast=True)
        future_reference_folder = dask_client.scatter(self.reference_folder, broadcast=True)
        future_configuration = dask_client.scatter(self.configuration, broadcast=True)

        for neoepitope in self.neoepitopes:
            logger.debug("Neoantigen: {}".format(neoepitope.to_json(indent=3)))
            futures.append(
                dask_client.submit(
                    NeoFoxEpitope.annotate_neoepitope,
                    neoepitope,
                    future_reference_folder,
                    future_configuration,
                    future_tcell_predictor,
                    future_self_similarity,
                    self.log_file_name,
                )
            )
        annotated_neoantigens = dask_client.gather(futures)
        end = time.time()
        logger.info(
            "Elapsed time for annotating {} neoepitopes {} seconds".format(
                len(self.neoepitopes), int(end - start)
            )
        )

        # close distributed resources
        del future_tcell_predictor
        del future_self_similarity
        del future_reference_folder
        del future_configuration

        return annotated_neoantigens

    @staticmethod
    def annotate_neoepitope(
        neoepitope: PredictedEpitope,
        reference_folder: ReferenceFolder,
        configuration: DependenciesConfiguration,
        tcell_predictor: TcellPrediction,
        self_similarity: SelfSimilarityCalculator,
        log_file_name: str,
    ):
        # the logs need to be initialised inside every dask job
        initialise_logs(log_file_name)
        logger.info("Starting neoepitope annotation with peptide={}".format(neoepitope.mutated_peptide))
        start = time.time()
        try:
            annotated_neoantigen = NeoepitopeAnnotator(
                reference_folder,
                configuration,
                tcell_predictor=tcell_predictor,
                self_similarity=self_similarity,
            ).get_annotated_neoepitope(neoepitope)
        except Exception as e:
            logger.error("Error processing neoantigen {}".format(neoepitope.to_dict()))
            raise e
        end = time.time()
        logger.info(
            "Elapsed time for annotating neoantigen for peptide={}: {} seconds".format(
                neoepitope.mutated_peptide, int(end - start))
        )
        return annotated_neoantigen

    def _conditional_expression_imputation(self) -> List[PredictedEpitope]:

        expression_annotator = ExpressionAnnotator()
        neoepitopes_transformed = []
        for neoepitope in self.neoepitopes:
            if neoepitope.patient_identifier is not None and neoepitope.patient_identifier != '':
                patient = self.patients[neoepitope.patient_identifier]
                neoepitope_transformed = neoepitope
                gene_expression = expression_annotator.get_gene_expression_annotation(
                    gene_name=neoepitope.gene, tcga_cohort=patient.tumor_type)
                neoepitope.imputed_gene_expression = gene_expression
                neoepitopes_transformed.append(neoepitope_transformed)
            else:
                neoepitopes_transformed.append(neoepitope)
        return neoepitopes_transformed


def initialise_logs(logfile, verbose=False):
    if logfile is not None:
        logzero.logfile(logfile)
    # TODO: this does not work
    if verbose:
        logzero.loglevel(logging.INFO)
    else:
        logzero.loglevel(logging.WARN)
