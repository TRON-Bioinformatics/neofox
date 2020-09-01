#!/usr/bin/env python
import logging
import sys
import os
import time
from typing import List

import logzero
from logzero import logger

from neofox import NEOFOX_LOG_FILE_ENV
from neofox.aa_index.aa_index import AminoacidIndex
from neofox.annotation_resources.gtex.gtex import GTEx
from neofox.annotation_resources.nmer_frequency.nmer_frequency import AminoacidFrequency, FourmerFrequency
from neofox.annotator import NeoantigenAnnotator
from neofox.exceptions import NeofoxConfigurationException
from neofox.helpers.available_alleles import AvailableAlleles
from neofox.helpers.runner import Runner
from neofox.annotation_resources.provean.provean import ProveanAnnotator
from neofox.model.conversion import ModelConverter
from neofox.model.neoantigen import NeoantigenAnnotations, Neoantigen, Patient
from neofox.predictors.MixMHCpred.mixmhc2pred import MixMhc2Pred
from neofox.predictors.MixMHCpred.mixmhcpred import MixMHCpred
from neofox.predictors.Tcell_predictor.tcellpredictor_wrapper import TcellPrediction
from neofox.predictors.dissimilarity_garnish.dissimilaritycalculator import DissimilarityCalculator
from neofox.predictors.neoag.neoag_gbm_model import NeoagCalculator
from neofox.predictors.neoantigen_fitness.neoantigen_fitness import NeoantigenFitnessCalculator
from neofox.predictors.netmhcpan4.combine_netmhcIIpan_pred_multiple_binders import BestAndMultipleBinderMhcII
from neofox.predictors.netmhcpan4.combine_netmhcpan_pred_multiple_binders import BestAndMultipleBinder
from neofox.references.references import ReferenceFolder, DependenciesConfiguration
from neofox.annotation_resources.uniprot.uniprot import Uniprot
from neofox.self_similarity.self_similarity import SelfSimilarityCalculator


class NeoFox:

    def __init__(self, neoantigens: List[Neoantigen], patient_id: str, patients: List[Patient]):

        logfile = os.environ.get(NEOFOX_LOG_FILE_ENV)
        if logfile is not None:
            logzero.logfile(logfile)
        # TODO: this does not work
        logzero.loglevel(logging.DEBUG)
        logger.info("Loading data...")
        if neoantigens is None or patients is None:
            raise NeofoxConfigurationException("Missing input data to run Neofox")
        self.neoantigens = neoantigens
        self.patients = {patient.identifier: patient for patient in patients}
        # TODO: avoid overriding patient id parameter
        for n in self.neoantigens:
            if n.patient_identifier is None:
                n.patient_identifier = patient_id

        self.references = ReferenceFolder()
        configuration = DependenciesConfiguration()
        runner = Runner()

        # resources with external dependencies (files or binaries)
        self.dissimilarity_calculator = DissimilarityCalculator(
            runner=runner, configuration=configuration, proteome_db=self.references.proteome_db)
        self.neoantigen_fitness_calculator = NeoantigenFitnessCalculator(
            runner=runner, configuration=configuration, iedb=self.references.iedb)
        self.neoag_calculator = NeoagCalculator(runner=runner, configuration=configuration)
        self.predII = BestAndMultipleBinderMhcII(runner=runner, configuration=configuration)
        self.predpresentation2 = MixMhc2Pred(runner=runner, configuration=configuration)
        self.pred = BestAndMultipleBinder(runner=runner, configuration=configuration)
        self.predpresentation = MixMHCpred(runner=runner, configuration=configuration)
        self.available_alleles = AvailableAlleles(self.references)
        self.provean_annotator = ProveanAnnotator(provean_file=self.references.prov_scores_mapped3)
        self.uniprot = Uniprot(self.references.uniprot)
        self.gtex = GTEx()
        self.aa_frequency = AminoacidFrequency()
        self.fourmer_frequency = FourmerFrequency()
        self.aa_index = AminoacidIndex()
        self.tcell_predictor = TcellPrediction()
        self.self_similarity = SelfSimilarityCalculator()
        logger.info("Data loaded")

    def get_annotations(self) -> List[NeoantigenAnnotations]:
        """
        Loads epitope data (if file has been not imported to R; colnames need to be changed), adds data to class that are needed to calculate,
        calls epitope class --> determination of epitope properties,
        write to txt file
        """
        logger.info("Starting NeoFox annotations...")
        epitope_annotator = NeoantigenAnnotator(
            provean_annotator=self.provean_annotator,
            uniprot=self.uniprot,
            dissimilarity_calculator=self.dissimilarity_calculator,
            neoantigen_fitness_calculator=self.neoantigen_fitness_calculator,
            neoag_calculator=self.neoag_calculator,
            netmhcpan2=self.predII,
            mixmhc2=self.predpresentation2,
            netmhcpan=self.pred,
            mixmhc=self.predpresentation,
            available_alleles=self.available_alleles,
            gtex=self.gtex,
            aa_frequency=self.aa_frequency,
            fourmer_frequency=self.fourmer_frequency,
            aa_index=self.aa_index,
            tcell_predictor=self.tcell_predictor,
            self_similarity=self.self_similarity
        )
        # feature calculation for each epitope
        annotations = []
        for neoantigen in self.neoantigens:
            logger.info("Annotating neoantigen...")
            start = time.time()
            patient = self.patients.get(neoantigen.patient_identifier)
            logger.debug("Neoantigen: {}".format(neoantigen.to_json(indent=3)))
            logger.debug("Patient: {}".format(patient.to_json(indent=3)))
            neoantigen_annotation = epitope_annotator.get_annotation(neoantigen, patient)
            end = time.time()
            logger.info("Elapsed time for annotating neoantigen {} seconds".format(int(end - start)))
            annotations.append(neoantigen_annotation)
        return annotations
