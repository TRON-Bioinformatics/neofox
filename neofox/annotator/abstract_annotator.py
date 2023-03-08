from abc import ABC

from neofox.MHC_predictors.netmhcpan.combine_netmhcpan_pred_multiple_binders import BestAndMultipleBinder
from neofox.annotation_resources.uniprot.uniprot import Uniprot
from neofox.helpers.blastp_runner import BlastpRunner
from neofox.helpers.epitope_helper import EpitopeHelper
from neofox.helpers.runner import Runner
from neofox.model.factories import AnnotationFactory
from neofox.model.neoantigen import PredictedEpitope, Neoantigen, Patient
from neofox.published_features.Tcell_predictor.tcellpredictor_wrapper import TcellPrediction
from neofox.published_features.differential_binding.amplitude import Amplitude
from neofox.published_features.differential_binding.differential_binding import DifferentialBinding
from neofox.published_features.dissimilarity_garnish.dissimilaritycalculator import DissimilarityCalculator
from neofox.published_features.hex.hex import Hex
from neofox.published_features.iedb_immunogenicity.iedb import IEDBimmunogenicity
from neofox.published_features.neoantigen_fitness.neoantigen_fitness import NeoantigenFitnessCalculator
from neofox.published_features.priority_score import PriorityScore
from neofox.published_features.self_similarity.self_similarity import SelfSimilarityCalculator
from neofox.references.references import ORGANISM_HOMO_SAPIENS, ReferenceFolder, DependenciesConfiguration


class AbstractAnnotator(ABC):

    def __init__(
            self,
            references: ReferenceFolder,
            configuration: DependenciesConfiguration,
            tcell_predictor: TcellPrediction,
            self_similarity: SelfSimilarityCalculator
    ):
        """class to annotate neoantigens"""
        self.runner = Runner()
        self.configuration = configuration
        self.tcell_predictor = tcell_predictor
        self.self_similarity = self_similarity
        self.organism = references.organism

        # NOTE: this one loads a big file, but it is faster loading it multiple times than passing it around
        self.uniprot = Uniprot(references.uniprot_pickle)

        # initialise proteome and IEDB BLASTP runners
        self.proteome_blastp_runner = BlastpRunner(
            runner=self.runner, configuration=configuration,
            database=references.get_proteome_database())
        self.iedb_blastp_runner = BlastpRunner(
            runner=self.runner, configuration=configuration,
            database=references.get_iedb_database())

        # NOTE: these resources do not read any file thus can be initialised fast
        self.dissimilarity_calculator = DissimilarityCalculator(proteome_blastp_runner=self.proteome_blastp_runner)
        self.neoantigen_fitness_calculator = NeoantigenFitnessCalculator(iedb_blastp_runner=self.iedb_blastp_runner)
        self.differential_binding = DifferentialBinding()
        self.priority_score_calculator = PriorityScore()
        self.iedb_immunogenicity = IEDBimmunogenicity()
        self.amplitude = Amplitude()
        self.hex = Hex(references=references)

    def get_additional_annotations_neoepitope_mhci(
            self, epitope: PredictedEpitope, neoantigen: Neoantigen = None) -> PredictedEpitope:

        if neoantigen is not None:
            gene = neoantigen.gene
            vaf_tumor_dna = neoantigen.dna_variant_allele_frequency
            vaf_tumor_rna = neoantigen.rna_variant_allele_frequency
            transcript_exp = neoantigen.rna_expression
        else:
            gene = epitope.gene
            vaf_tumor_dna = epitope.dna_variant_allele_frequency
            vaf_tumor_rna = epitope.rna_variant_allele_frequency
            transcript_exp = epitope.rna_expression

        epitope.neofox_annotations.annotations.extend(
            BestAndMultipleBinder.get_annotations_epitope_mhci(epitope=epitope) +
            self.amplitude.get_annotations_epitope_mhci(epitope=epitope)
        )

        # NOTE: this extend() call cannot be joined with the previous as some of the previous annotations are expected
        epitope.neofox_annotations.annotations.extend(
            self.neoantigen_fitness_calculator.get_annotations_epitope_mhci(epitope=epitope) +
            self.differential_binding.get_annotations_epitope_mhci(epitope=epitope) +
            self.self_similarity.get_annotations_epitope_mhci(epitope=epitope) +
            self.uniprot.get_annotations_epitope(epitope=epitope) +
            self.dissimilarity_calculator.get_annotations_epitope(epitope=epitope)
        )

        epitope.neofox_annotations.annotations.extend(self.tcell_predictor.get_annotations_epitope_mhci(
            epitope=epitope, gene=gene))

        num_mismatches = EpitopeHelper.number_of_mismatches(
            epitope_wild_type=epitope.wild_type_peptide, epitope_mutation=epitope.mutated_peptide, )
        epitope.neofox_annotations.annotations.append(AnnotationFactory.build_annotation(
            value=num_mismatches,
            name='number_of_mismatches'))

        epitope.neofox_annotations.annotations.extend(
            self.priority_score_calculator.get_annotations_epitope_mhci(
                epitope=epitope, vaf_rna=vaf_tumor_rna, vaf_tumor=vaf_tumor_dna, transcript_exp=transcript_exp))

        if self.organism == ORGANISM_HOMO_SAPIENS:
            epitope.neofox_annotations.annotations.extend(
                self.iedb_immunogenicity.get_annotations_epitope_mhci(epitope=epitope) +
                self.hex.get_annotations_epitope(epitope=epitope))

        return epitope

    def get_additional_annotations_neoepitope_mhcii(self, epitope: PredictedEpitope) -> PredictedEpitope:

        epitope.neofox_annotations.annotations.extend(
            self.amplitude.get_annotations_epitope_mhcii(epitope=epitope) +
            self.neoantigen_fitness_calculator.get_annotations_epitope_mhcii(epitope=epitope) +
            self.self_similarity.get_annotations_epitope_mhcii(epitope=epitope) +
            self.uniprot.get_annotations_epitope(epitope=epitope) +
            self.dissimilarity_calculator.get_annotations_epitope(epitope=epitope))

        if self.organism == ORGANISM_HOMO_SAPIENS:

            epitope.neofox_annotations.annotations.extend(
                self.iedb_immunogenicity.get_annotations_epitope_mhcii(epitope=epitope) +
                self.hex.get_annotations_epitope(epitope=epitope))

        return epitope