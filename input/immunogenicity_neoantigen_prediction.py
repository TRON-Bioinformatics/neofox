#!/usr/bin/env python

from logzero import logger

from input.aa_index.aa_index import AaIndex
from input.IEDB_Immunogenicity.predict_immunogenicity_simple import IEDBimmunogenicity
from input.annotation_resources.nmer_frequency.nmer_frequency import AminoacidFrequency, FourmerFrequency
from input.epitope_annotator import EpitopeAnnotator
from input.annotation_resources.gtex.gtex import GTEx
from input.helpers import data_import
from input.helpers.available_alleles import AvailableAlleles
from input.helpers.properties_manager import PATIENT_ID
from input.helpers.runner import Runner
from input.literature_features.differential_binding import DifferentialBinding
from input.literature_features.expression import Expression
from input.literature_features.priority_score import PriorityScore
from input.new_features.conservation_scores import ProveanAnnotator
from input.predictors.MixMHCpred.mixmhc2pred import MixMhc2Pred
from input.predictors.MixMHCpred.mixmhcpred import MixMHCpred
from input.predictors.Tcell_predictor.tcellpredictor_wrapper import TcellPrediction
from input.predictors.dissimilarity_garnish.dissimilaritycalculator import DissimilarityCalculator
from input.predictors.neoag.neoag_gbm_model import NeoagCalculator
from input.predictors.neoantigen_fitness.neoantigen_fitness import NeoantigenFitnessCalculator
from input.predictors.netmhcpan4.combine_netmhcIIpan_pred_multiple_binders import BestAndMultipleBinderMhcII
from input.predictors.netmhcpan4.combine_netmhcpan_pred_multiple_binders import BestAndMultipleBinder
from input.references import ReferenceFolder, DependenciesConfiguration
from input.annotation_resources.uniprot.uniprot import Uniprot


class ImmunogenicityNeoantigenPredictionToolbox:

    def __init__(self, icam_file, patient_id, patients_file):

        self.patient_id = patient_id

        self.references = ReferenceFolder()
        configuration = DependenciesConfiguration()
        runner = Runner()

        # resources without external dependencies
        self.gtex = GTEx()
        self.uniprot = Uniprot(self.references.uniprot)
        self.aa_frequency = AminoacidFrequency()
        self.fourmer_frequency = FourmerFrequency()
        self.aa_index = AaIndex()
        self.differential_binding = DifferentialBinding()
        self.expression_calculator = Expression()
        self.priority_score_calcualtor = PriorityScore()
        self.tcell_predictor = TcellPrediction()

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
        self.iedb_immunogenicity = IEDBimmunogenicity()
        self.available_alleles = AvailableAlleles(self.references)

        # import epitope data
        self.header, self.rows = data_import.import_dat_icam(icam_file)
        self.patients = {patient.identifier: patient for patient in data_import.import_patients_data(patients_file)}

        # TODO: remove once we are loading data into models
        if "+-13_AA_(SNV)_/_-15_AA_to_STOP_(INDEL)" in self.header:
            self.header, self.rows = data_import.change_col_names(header=self.header, data=self.rows)

        # adds patient to the table
        self.header.append(PATIENT_ID)
        logger.debug(self.patient_id)

        # TODO: when we are moving away from the icam table we will have a patient id for each neoantigen and
        # TODO: we will be able to pass the whole list of patients below
        self.tissue = self.patients.get(self.patient_id).tissue
        for row in self.rows:
            row.append(str(self.patient_id))

        self.provean_annotator = ProveanAnnotator(provean_file=self.references.prov_scores_mapped3,
                                                  header_epitopes=self.header, epitopes=self.rows)

    def get_annotations(self):
        """
        Loads epitope data (if file has been not imported to R; colnames need to be changed), adds data to class that are needed to calculate,
        calls epitope class --> determination of epitope properties,
        write to txt file
        """
        epitope_annotator = EpitopeAnnotator(
            provean_annotator=self.provean_annotator,
            gtex=self.gtex,
            uniprot=self.uniprot,
            aa_frequency=self.aa_frequency,
            fourmer_frequency=self.fourmer_frequency,
            aa_index=self.aa_index,
            dissimilarity_calculator=self.dissimilarity_calculator,
            neoantigen_fitness_calculator=self.neoantigen_fitness_calculator,
            neoag_calculator=self.neoag_calculator,
            predII=self.predII,
            predpresentation2=self.predpresentation2,
            pred=self.pred,
            predpresentation=self.predpresentation,
            tcell_predictor=self.tcell_predictor,
            iedb_immunogenicity=self.iedb_immunogenicity,
            differential_binding=self.differential_binding,
            expression_calculator=self.expression_calculator,
            priority_score_calculator=self.priority_score_calcualtor,
            available_alleles=self.available_alleles,
            patients=self.patients)
        # feature calculation for each epitope
        annotations = []
        for row in self.rows:
            # TODO: move this initialisation out of the loop once the properties have been refactored out
            annotation = epitope_annotator.get_annotation(self.header, row, self.patient_id, self.tissue)
            annotations.append(annotation)
        return annotations, self.header
