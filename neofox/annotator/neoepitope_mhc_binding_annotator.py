from typing import Tuple

from neofox.MHC_predictors.MixMHCpred.mixmhc2pred import MixMHC2pred
from neofox.MHC_predictors.MixMHCpred.mixmhcpred import MixMHCpred
from neofox.MHC_predictors.netmhcpan.combine_netmhcIIpan_pred_multiple_binders import BestAndMultipleBinderMhcII
from neofox.MHC_predictors.netmhcpan.combine_netmhcpan_pred_multiple_binders import BestAndMultipleBinder
from neofox.MHC_predictors.netmhcpan.netmhcIIpan_prediction import NetMhcIIPanPredictor
from neofox.MHC_predictors.netmhcpan.netmhcpan_prediction import NetMhcPanPredictor
from neofox.MHC_predictors.prime import Prime
from neofox.annotation_resources.uniprot.uniprot import Uniprot
from neofox.helpers.blastp_runner import BlastpRunner
from neofox.helpers.epitope_helper import EpitopeHelper
from neofox.helpers.runner import Runner
from neofox.model.factories import AnnotationFactory
from neofox.model.mhc_parser import MhcParser
from neofox.model.neoantigen import Neoantigen, Patient, PredictedEpitope
from neofox.references.references import DependenciesConfiguration, AvailableAlleles, ReferenceFolder, \
    ORGANISM_HOMO_SAPIENS


class NeoepitopeMhcBindingAnnotator:

    def __init__(self, references: ReferenceFolder, configuration: DependenciesConfiguration,
                 uniprot: Uniprot, proteome_blastp_runner: BlastpRunner):
        """class to annotate neoantigens"""
        self.runner = Runner()
        self.configuration = configuration
        self.proteome_db = references.proteome_db
        self.available_alleles = references.get_available_alleles()
        self.organism = references.organism
        self.uniprot = uniprot
        self.proteome_blastp_runner = proteome_blastp_runner

        self.mhc_database = references.get_mhc_database()
        self.mhc_parser = MhcParser.get_mhc_parser(self.mhc_database)
        self.netmhcpan = NetMhcPanPredictor(
            runner=self.runner, configuration=configuration, mhc_parser=self.mhc_parser,
            blastp_runner=self.proteome_blastp_runner)
        self.netmhc2pan = NetMhcIIPanPredictor(
            runner=self.runner, configuration=configuration, mhc_parser=self.mhc_parser,
            blastp_runner=self.proteome_blastp_runner)
        self.mixmhcpred = MixMHCpred(self.runner, self.configuration, self.mhc_parser)
        self.mixmhc2pred = MixMHC2pred(self.runner, self.configuration, self.mhc_parser)
        self.prime = Prime(self.runner, self.configuration, self.mhc_parser)

    def get_mhc_binding_annotations(self, neoepitope: PredictedEpitope) -> PredictedEpitope:

        has_mhc1 = neoepitope.allele_mhc_i is not None and neoepitope.allele_mhc_i.name != ''
        has_mhc2 = neoepitope.isoform_mhc_i_i is not None and neoepitope.isoform_mhc_i_i.name != ''

        if has_mhc1:
            # MHC I epitope
            annotated_neoepitope = self._run_netmhcpan(neoepitope=neoepitope)
            if annotated_neoepitope:
                if self.configuration.mix_mhc_pred and self.organism == ORGANISM_HOMO_SAPIENS:
                    mixmhcpred_neoepitope, mixmhcpred_neoepitope_wt = self._run_mixmhcpred(neoepitope=neoepitope)
                    if mixmhcpred_neoepitope:
                        annotated_neoepitope = AnnotationFactory.annotate_epitope(
                            epitope=annotated_neoepitope,
                            paired_epitope=mixmhcpred_neoepitope,
                            annotation_name=MixMHCpred.ANNOTATION_PREFIX)
                        annotated_neoepitope = AnnotationFactory.annotate_epitope(
                            epitope=annotated_neoepitope,
                            paired_epitope=mixmhcpred_neoepitope_wt,
                            annotation_name=MixMHCpred.ANNOTATION_PREFIX_WT)
                    if self.configuration.prime:
                        prime_neoepitope, prime_neoepitope_wt = self._run_prime(neoepitope=neoepitope)
                        if prime_neoepitope:
                            annotated_neoepitope = AnnotationFactory.annotate_epitope(
                                epitope=annotated_neoepitope,
                                paired_epitope=prime_neoepitope,
                                annotation_name=Prime.ANNOTATION_PREFIX)
                            annotated_neoepitope = AnnotationFactory.annotate_epitope(
                                epitope=annotated_neoepitope,
                                paired_epitope=prime_neoepitope_wt,
                                annotation_name=Prime.ANNOTATION_PREFIX_WT)
        elif has_mhc2:
            # MHC II epitope
            annotated_neoepitope = self._run_netmhc2pan(neoepitope=neoepitope)
            if annotated_neoepitope:
                if self.configuration.mix_mhc2_pred and self.organism == ORGANISM_HOMO_SAPIENS:
                    mixmhc2pred_neoepitope, mixmhc2pred_neoepitope_wt = self._run_mixmhc2pred(neoepitope=neoepitope)
                    if mixmhc2pred_neoepitope:
                        annotated_neoepitope = AnnotationFactory.annotate_epitope(
                            epitope=annotated_neoepitope,
                            paired_epitope=mixmhc2pred_neoepitope,
                            annotation_name=MixMHC2pred.ANNOTATION_PREFIX)
                        annotated_neoepitope = AnnotationFactory.annotate_epitope(
                            epitope=annotated_neoepitope,
                            paired_epitope=mixmhc2pred_neoepitope_wt,
                            annotation_name=MixMHC2pred.ANNOTATION_PREFIX_WT)
        else:
            raise ValueError("Neoepitope without neither MHC I allele or MHC II isoform")

        return annotated_neoepitope

    def _run_netmhcpan(self, neoepitope: PredictedEpitope) -> PredictedEpitope:
        # runs NetMHCpan in peptide mode over the mutated and WT separately and merges it back in one
        # predicted epitope
        annotated_neoepitope = neoepitope

        netmhcpan_allele = self.mhc_parser.get_netmhcpan_representation(neoepitope.allele_mhc_i)
        if netmhcpan_allele in self.available_alleles.get_available_mhc_i():
            mutated_epitope = self.netmhcpan.mhc_prediction_peptide(
                sequence=neoepitope.mutated_peptide, alleles=netmhcpan_allele)
            annotated_neoepitope.affinity_mutated = mutated_epitope.affinity_mutated
            annotated_neoepitope.rank_mutated = mutated_epitope.rank_mutated
            wt_epitope = self.netmhcpan.mhc_prediction_peptide(
                sequence=neoepitope.wild_type_peptide, alleles=netmhcpan_allele)
            annotated_neoepitope.affinity_wild_type = wt_epitope.affinity_mutated
            annotated_neoepitope.rank_wild_type = wt_epitope.rank_mutated
        return annotated_neoepitope

    def _run_netmhc2pan(self, neoepitope: PredictedEpitope) -> PredictedEpitope:
        annotated_neoepitope = neoepitope

        netmhc2pan_allele = self.mhc_parser.get_netmhc2pan_representation(neoepitope.isoform_mhc_i_i)
        if netmhc2pan_allele in self.available_alleles.get_available_mhc_ii():
            mutated_epitope = self.netmhc2pan.mhc2_prediction_peptide(
                sequence=neoepitope.mutated_peptide,
                mhc2_isoform=neoepitope.isoform_mhc_i_i)
            annotated_neoepitope.affinity_mutated = mutated_epitope.affinity_mutated
            annotated_neoepitope.rank_mutated = mutated_epitope.rank_mutated
            wt_epitope = self.netmhc2pan.mhc2_prediction_peptide(
                sequence=neoepitope.wild_type_peptide,
                mhc2_isoform=neoepitope.isoform_mhc_i_i)
            annotated_neoepitope.affinity_wild_type = wt_epitope.affinity_mutated
            annotated_neoepitope.rank_wild_type = wt_epitope.rank_mutated
        return annotated_neoepitope

    def _run_mixmhcpred(self, neoepitope: PredictedEpitope) -> Tuple[PredictedEpitope, PredictedEpitope]:
        mutated_epitope = self.mixmhcpred.run_peptide(
            peptide=neoepitope.mutated_peptide, allele=neoepitope.allele_mhc_i)
        wt_epitope = self.mixmhcpred.run_peptide(
            peptide=neoepitope.wild_type_peptide, allele=neoepitope.allele_mhc_i)
        return mutated_epitope, wt_epitope

    def _run_prime(self, neoepitope: PredictedEpitope) -> Tuple[PredictedEpitope, PredictedEpitope]:
        mutated_epitope = self.prime.run_peptide(peptide=neoepitope.mutated_peptide, allele=neoepitope.allele_mhc_i)
        wt_epitope = self.prime.run_peptide(
            peptide=neoepitope.wild_type_peptide, allele=neoepitope.allele_mhc_i)
        return mutated_epitope, wt_epitope

    def _run_mixmhc2pred(self, neoepitope: PredictedEpitope) -> Tuple[PredictedEpitope, PredictedEpitope]:
        mutated_epitope = self.mixmhc2pred.run_peptide(
            isoform=neoepitope.isoform_mhc_i_i, peptide=neoepitope.mutated_peptide)
        wt_epitope = self.mixmhc2pred.run_peptide(
            peptide=neoepitope.wild_type_peptide, isoform=neoepitope.isoform_mhc_i_i)
        return mutated_epitope, wt_epitope
