#!/usr/bin/env python
import logging
import sys
import os
import time
from typing import List

import logzero
from logzero import logger

from neofox import NEOFOX_LOG_FILE_ENV
from neofox.annotator import NeoantigenAnnotator
from neofox.helpers.available_alleles import AvailableAlleles
from neofox.helpers.runner import Runner
from neofox.annotation_resources.provean.provean import ProveanAnnotator
from neofox.model.conversion import ModelConverter
from neofox.model.neoantigen import NeoantigenAnnotations
from neofox.predictors.MixMHCpred.mixmhc2pred import MixMhc2Pred
from neofox.predictors.MixMHCpred.mixmhcpred import MixMHCpred
from neofox.predictors.dissimilarity_garnish.dissimilaritycalculator import DissimilarityCalculator
from neofox.predictors.neoag.neoag_gbm_model import NeoagCalculator
from neofox.predictors.neoantigen_fitness.neoantigen_fitness import NeoantigenFitnessCalculator
from neofox.predictors.netmhcpan4.combine_netmhcIIpan_pred_multiple_binders import BestAndMultipleBinderMhcII
from neofox.predictors.netmhcpan4.combine_netmhcpan_pred_multiple_binders import BestAndMultipleBinder
from neofox.references.references import ReferenceFolder, DependenciesConfiguration
from neofox.annotation_resources.uniprot.uniprot import Uniprot


class NeoFox:

    def __init__(self, icam_file, patient_id, patients_file):

        logfile = os.environ.get(NEOFOX_LOG_FILE_ENV)
        if logfile is not None:
            logzero.logfile(logfile)
        # TODO: this does not work
        logzero.loglevel(logging.DEBUG)
        logger.info("Loading data...")
        self.patient_id = patient_id
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

        # import epitope data
        self.neoantigens = ModelConverter.parse_icam_file(icam_file)
        self.patients = {patient.identifier: patient for patient in ModelConverter.parse_patients_file(patients_file)}

        # TODO: avoid overriding patient id parameter
        for n in self.neoantigens:
            if n.patient_identifier is None:
                n.patient_identifier = self.patient_id
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
            available_alleles=self.available_alleles)
        # feature calculation for each epitope
        annotations = []
        for neoantigen in self.neoantigens:
            logger.info("Annotating neoantigen...")
            start = time.time()
            patient = self.patients.get(neoantigen.patient_identifier)
            logger.debug("Neoantigen: {}".format(neoantigen.to_json(indent=3)))
            logger.debug("Patient: {}".format(patient.to_json(indent=3)))
            annotation = epitope_annotator.get_annotation(neoantigen, patient)
            end = time.time()
            logger.info("Elapsed time for annotating neoantigen {} seconds".format(int(end - start)))
            annotations.append(annotation)
        return annotations

    @staticmethod
    def write_to_file_sorted(annotations: NeoantigenAnnotations, output_file=None):
        """Transforms dictionary (property --> epitopes). To one unit (epitope) corresponding values are concentrated in one list
        and printed ';' separated."""
        logger.info("Writing results...")
        if output_file is not None:
            sys.stdout = open(output_file, 'w')     # redirects stdout if file provided
        transformed_annotations = {}
        for neoantigen in annotations:
            for key in neoantigen:
                if key not in transformed_annotations:
                    # keys are are feautres; values: list of feature values associated with mutated peptide sequence
                    transformed_annotations[key] = [neoantigen[key]]
                else:
                    transformed_annotations[key].append(neoantigen[key])

        features_names = []
        for key in transformed_annotations:
            if key not in header:
                features_names.append(key)
        features_names.sort()
        header.extend(features_names)
        print("\t".join(header))
        for i in range(len(transformed_annotations["patient_identifier"])):  # NOTE: this has nothing to do with "patient_identifier" field
            z = [NeoFox.fetch_annotation(transformed_annotations, col, i) for col in header]
            print("\t".join(z))
        sys.stdout.flush()  # this is required to avoid half written files
        logger.info("Finished writing")

    @staticmethod
    def fetch_annotation(annotations, column, index):
        # TODO: this is a compromise solution that will go away soon
        result = "NA"
        try:
            result = str(annotations[column][index])
        except IndexError:
            pass
        return result
