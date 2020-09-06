#!/usr/bin/env python
from typing import List

import neofox.predictors.netmhcpan4.multiple_binders as multiple_binders
import neofox.predictors.netmhcpan4.netmhcpan_prediction as netmhcpan_prediction
from neofox.helpers import intermediate_files
from neofox.helpers.epitope_helper import EpitopeHelper
from neofox.model.neoantigen import Annotation
from neofox.model.wrappers import AnnotationFactory
import neofox.helpers.casting as casting


class BestAndMultipleBinder:
    def __init__(self, runner, configuration):
        """
        :type runner: neofox.helpers.runner.Runner
        :type configuration: neofox.references.DependenciesConfiguration
        """
        self.runner = runner
        self.configuration = configuration
        self.mean_type = ["arithmetic", "harmonic", "geometric"]
        self._initialise()

    def _initialise(self):
        self.MHC_score_all_epitopes = []
        self.MHC_score_top10 = []
        self.MHC_score_best_per_alelle = []
        self.MHC_number_strong_binders = None
        self.MHC_number_weak_binders = None
        self.MHC_epitope_seqs = None
        self.MHC_epitope_scores = None
        self.MHC_epitope_alleles = None
        self.best4_mhc_score = None
        self.best4_mhc_epitope = None
        self.best4_mhc_allele = None
        self.best4_mhc_position = None
        self.directed_to_TCR = None
        self.best4_affinity = None
        self.best4_affinity_epitope = "-"
        self.best4_affinity_allele = None
        self.best4_affinity_position = None
        self.best4_affinity_directed_to_TCR = None
        self.epitope_affinities = None
        self.generator_rate = None
        self.mhcI_score_9mer = None
        self.mhcI_score_allele_9mer = None
        self.mhcI_score_position_9mer = None
        self.mhcI_score_epitope_9mer = "-"
        self.mhcI_affinity_9mer = None
        self.mhcI_affinity_allele_9mer = None
        self.mhcI_affinity_position_9mer = None
        self.mhcI_affinity_epitope_9mer = "-"
        # WT features
        self.MHC_epitope_scores_WT = None
        self.MHC_epitope_seqs_WT = None
        self.MHC_epitope_alleles_WT = None
        self.MHC_score_top10_WT = []
        self.MHC_score_all_epitopes_WT = []
        self.MHC_score_best_per_alelle_WT = []
        self.MHC_number_strong_binders_WT = None
        self.MHC_number_weak_binders_WT = None
        self.best4_mhc_score_WT = None
        self.best4_mhc_epitope_WT = None
        self.best4_mhc_allele_WT = None
        self.best4_affinity_WT = None
        self.best4_affinity_epitope_WT = None
        self.best4_affinity_allele_WT = None
        self.epitope_affinities_WT = None
        self.generator_rate_WT = None
        self.mhcI_score_9mer_WT = None
        self.mhcI_score_allele_9mer_WT = None
        self.mhcI_score_epitope_9mer_WT = None
        self.mhcI_affinity_9mer_WT = None
        self.mhcI_affinity_allele_9mer_WT = None
        self.mhcI_affinity_epitope_9mer_WT = None

    def mhc_mb_score_best_per_allele(self, tuple_best_per_allele):
        '''returns list of multiple binding scores for mhcII considering best epitope per allele, applying different types of means (harmonic ==> PHRB-II, Marty et al).
        2 copies of DRA - DRB1 --> consider this gene 2x when averaging mhcii binding scores
        '''
        multbind = multiple_binders.MultipleBinding()
        tuple_best_per_allele_new = list(tuple_best_per_allele)
        mean_best_per_per_allele = [None, None, None]
        if len(tuple_best_per_allele_new) == 6:
            mean_best_per_per_allele = multbind.get_means(tuple_best_per_allele_new)
        return mean_best_per_per_allele

    def run(self, xmer_wt, xmer_mut, alleles, set_available_mhc):
        """
        predicts MHC epitopes; returns on one hand best binder and on the other hand multiple binder analysis is performed
        """
        ### PREDICTION FOR MUTATED SEQUENCE
        self._initialise()
        tmp_prediction = intermediate_files.create_temp_file(prefix="netmhcpanpred_", suffix=".csv")
        np = netmhcpan_prediction.NetMhcPanPredictor(runner=self.runner, configuration=self.configuration)
        mb = multiple_binders.MultipleBinding()
        tmp_fasta = intermediate_files.create_temp_fasta(sequences=[xmer_mut], prefix="tmp_singleseq_")
        # print alleles
        np.mhc_prediction(alleles, set_available_mhc, tmp_fasta, tmp_prediction)

        position_xmer = np.mut_position_xmer_seq(xmer_mut=xmer_mut, xmer_wt=xmer_wt)
        preds = np.filter_binding_predictions(position_xmer=position_xmer, tmppred=tmp_prediction)

        # multiple binding
        list_tups = mb.generate_epi_tuple(preds)
        self.MHC_epitope_scores = "/".join([str(tup[0]) for tup in list_tups])
        self.epitope_affinities = "/".join([str(tup[1]) for tup in list_tups])
        self.MHC_epitope_seqs = "/".join([tup[2] for tup in list_tups])
        self.MHC_epitope_alleles = "/".join([tup[3] for tup in list_tups])
        top10 = mb.extract_top10_epis(list_tups)
        best_per_alelle = mb.extract_best_epi_per_alelle(list_tups, alleles)
        all = mb.scores_to_list(list_tups)
        all_affinities = mb.affinities_to_list(list_tups)
        top10 = mb.scores_to_list(top10)
        self.MHC_score_top10 = mb.get_means(top10)
        best_per_alelle = mb.scores_to_list(best_per_alelle)
        self.MHC_score_all_epitopes = mb.get_means(all)
        self.MHC_score_best_per_alelle = self.mhc_mb_score_best_per_allele(best_per_alelle)
        self.MHC_number_strong_binders = mb.determine_number_of_binders(all, 1)
        self.MHC_number_weak_binders = mb.determine_number_of_binders(all, 2)
        # best prediction
        best_epi = np.minimal_binding_score(preds)
        self.best4_mhc_score = casting.to_float(np.add_best_epitope_info(best_epi, "%Rank"))
        self.best4_mhc_epitope = np.add_best_epitope_info(best_epi, "Peptide")
        self.best4_mhc_allele = np.add_best_epitope_info(best_epi, "HLA")
        self.best4_mhc_position = np.add_best_epitope_info(best_epi, "Pos")
        self.directed_to_TCR = np.mutation_in_loop(position_xmer_list=position_xmer, epitope_tuple=best_epi)
        best_epi_affinity = np.minimal_binding_score(preds, rank=False)
        self.best4_affinity = casting.to_float(np.add_best_epitope_info(best_epi_affinity, "Aff(nM)"))
        self.best4_affinity_epitope = np.add_best_epitope_info(best_epi_affinity, "Peptide")
        self.best4_affinity_allele = np.add_best_epitope_info(best_epi_affinity, "HLA")
        self.best4_affinity_position = np.add_best_epitope_info(best_epi_affinity, "Pos")
        self.best4_affinity_directed_to_TCR = np.mutation_in_loop(
            position_xmer_list=position_xmer, epitope_tuple=best_epi_affinity)
        # multiple binding based on affinity
        self.generator_rate = mb.determine_number_of_binders(list_scores=all_affinities, threshold=50)
        # best predicted epitope of length 9
        preds_9mer = np.filter_for_9mers(preds)
        best_9mer = np.minimal_binding_score(preds_9mer)
        best_9mer_affinity = np.minimal_binding_score(preds_9mer, rank=False)
        self.mhcI_score_9mer = np.add_best_epitope_info(best_9mer, "%Rank")
        self.mhcI_score_allele_9mer = np.add_best_epitope_info(best_9mer, "HLA")
        self.mhcI_score_position_9mer = np.add_best_epitope_info(best_9mer, "Pos")
        self.mhcI_score_epitope_9mer = np.add_best_epitope_info(best_9mer, "Peptide")
        self.mhcI_affinity_9mer = casting.to_float(np.add_best_epitope_info(best_9mer_affinity, "Aff(nM)"))
        self.mhcI_affinity_allele_9mer = np.add_best_epitope_info(best_9mer_affinity, "HLA")
        self.mhcI_affinity_position_9mer = np.add_best_epitope_info(best_9mer_affinity, "Pos")
        self.mhcI_affinity_epitope_9mer = np.add_best_epitope_info(best_9mer_affinity, "Peptide")

        ### PREDICTION FOR WT SEQUENCE
        tmp_prediction = intermediate_files.create_temp_file(prefix="netmhcpanpred_", suffix=".csv")
        np = netmhcpan_prediction.NetMhcPanPredictor(runner=self.runner, configuration=self.configuration)
        mb = multiple_binders.MultipleBinding()
        tmp_fasta = intermediate_files.create_temp_fasta(sequences=[xmer_wt],
                                                         prefix="tmp_singleseq_")
        np.mhc_prediction(alleles, set_available_mhc, tmp_fasta, tmp_prediction)
        preds = np.filter_binding_predictions(position_xmer=position_xmer, tmppred=tmp_prediction)
        # multiple binding
        list_tups = mb.generate_epi_tuple(preds)
        self.MHC_epitope_scores_WT = "/".join([str(tup[0]) for tup in list_tups])
        self.epitope_affinities_WT = "/".join([str(tup[1]) for tup in list_tups])
        self.MHC_epitope_seqs_WT = "/".join([tup[2] for tup in list_tups])
        self.MHC_epitope_alleles_WT = "/".join([tup[3] for tup in list_tups])
        top10 = mb.extract_top10_epis(list_tups)
        best_per_alelle = mb.extract_best_epi_per_alelle(list_tups, alleles)
        all = mb.scores_to_list(list_tups)
        all_affinities = mb.affinities_to_list(list_tups)
        top10 = mb.scores_to_list(top10)
        self.MHC_score_top10_WT = mb.get_means(top10)
        best_per_alelle = mb.scores_to_list(best_per_alelle)
        self.MHC_score_all_epitopes_WT = mb.get_means(all)
        self.MHC_score_best_per_alelle_WT = self.mhc_mb_score_best_per_allele(best_per_alelle)
        self.MHC_number_strong_binders_WT = mb.determine_number_of_binders(all, 1)
        self.MHC_number_weak_binders_WT = mb.determine_number_of_binders(all, 2)
        # best prediction
        best_epi = np.filter_for_WT_epitope_position(preds, self.best4_mhc_epitope,
                                            position_epi_xmer=self.best4_mhc_position)
        self.best4_mhc_score_WT = casting.to_float(np.add_best_epitope_info(best_epi, "%Rank"))
        self.best4_mhc_epitope_WT = np.add_best_epitope_info(best_epi, "Peptide")
        self.best4_mhc_allele_WT = np.add_best_epitope_info(best_epi, "HLA")

        best_epi_affinity = np.filter_for_WT_epitope_position(preds, self.best4_affinity_epitope,
                                                     position_epi_xmer=self.best4_affinity_position)
        self.best4_affinity_WT = casting.to_float(np.add_best_epitope_info(best_epi_affinity, "Aff(nM)"))
        self.best4_affinity_epitope_WT = np.add_best_epitope_info(best_epi_affinity, "Peptide")
        self.best4_affinity_allele_WT = np.add_best_epitope_info(best_epi_affinity, "HLA")
        self.generator_rate_WT = mb.determine_number_of_binders(list_scores=all_affinities, threshold=50)
        # best predicted epitope of length 9
        preds_9mer = np.filter_for_9mers(preds)
        best_9mer = np.filter_for_WT_epitope_position(preds_9mer, self.mhcI_score_epitope_9mer,
                                             position_epi_xmer=self.mhcI_score_position_9mer)
        best_9mer_affinity = np.filter_for_WT_epitope_position(preds_9mer, mut_seq=self.mhcI_affinity_epitope_9mer,
                                                      position_epi_xmer=self.mhcI_affinity_position_9mer)
        self.mhcI_score_9mer_WT = np.add_best_epitope_info(best_9mer, "%Rank")
        self.mhcI_score_allele_9mer_WT = np.add_best_epitope_info(best_9mer, "HLA")
        self.mhcI_score_epitope_9mer_WT = np.add_best_epitope_info(best_9mer, "Peptide")
        self.mhcI_affinity_9mer_WT = casting.to_float(np.add_best_epitope_info(best_9mer_affinity, "Aff(nM)"))
        self.mhcI_affinity_allele_9mer_WT = np.add_best_epitope_info(best_9mer_affinity, "HLA")
        self.mhcI_affinity_epitope_9mer_WT = np.add_best_epitope_info(best_9mer_affinity, "Peptide")

    def get_annotations(self) -> List[Annotation]:
        annotations = [AnnotationFactory.build_annotation(value=self.best4_mhc_score, name="best%Rank_netmhcpan4"),
                       AnnotationFactory.build_annotation(value=self.best4_mhc_epitope, name="best_epitope_netmhcpan4"),
                       AnnotationFactory.build_annotation(value=self.best4_mhc_allele, name="bestHLA_allele_netmhcpan4"),
                       AnnotationFactory.build_annotation(value=self.directed_to_TCR, name="directed_to_TCR"),
                       AnnotationFactory.build_annotation(value=self.best4_affinity, name="best_affinity_netmhcpan4"),
                       AnnotationFactory.build_annotation(value=self.best4_affinity_epitope, name="best_affinity_epitope_netmhcpan4"),
                       AnnotationFactory.build_annotation(value=self.best4_affinity_allele, name="bestHLA_allele_affinity_netmhcpan4"),
                       AnnotationFactory.build_annotation(value=self.best4_affinity_directed_to_TCR, name="affinity_directed_to_TCR"),
                       AnnotationFactory.build_annotation(value=self.mhcI_score_9mer, name="best%Rank_netmhcpan4_9mer"),
                       AnnotationFactory.build_annotation(value=self.mhcI_score_epitope_9mer, name="best_epitope_netmhcpan4_9mer"),
                       AnnotationFactory.build_annotation(value=self.mhcI_score_allele_9mer, name="bestHLA_allele_netmhcpan4_9mer"),
                       AnnotationFactory.build_annotation(value=self.mhcI_affinity_9mer, name="best_affinity_netmhcpan4_9mer"),
                       AnnotationFactory.build_annotation(value=self.mhcI_affinity_allele_9mer, name="bestHLA_allele_affinity_netmhcpan4_9mer"),
                       AnnotationFactory.build_annotation(value=self.mhcI_affinity_epitope_9mer, name="best_affinity_epitope_netmhcpan4_9mer"),
                       AnnotationFactory.build_annotation(value=self.best4_affinity_WT, name="best_affinity_netmhcpan4_WT"),
                       AnnotationFactory.build_annotation(value=self.best4_affinity_epitope_WT, name="best_affinity_epitope_netmhcpan4_WT"),
                       AnnotationFactory.build_annotation(value=self.best4_affinity_allele_WT, name="bestHLA_allele_affinity_netmhcpan4_WT"),
                       AnnotationFactory.build_annotation(value=self.best4_mhc_score_WT, name="best%Rank_netmhcpan4_WT"),
                       AnnotationFactory.build_annotation(value=self.best4_mhc_epitope_WT, name="best_epitope_netmhcpan4_WT"),
                       AnnotationFactory.build_annotation(value=self.best4_mhc_allele_WT, name="bestHLA_allele_netmhcpan4_WT"),
                       AnnotationFactory.build_annotation(value=self.mhcI_score_9mer_WT, name="best%Rank_netmhcpan4_9mer_WT"),
                       AnnotationFactory.build_annotation(value=self.mhcI_score_epitope_9mer_WT, name="best_epitope_netmhcpan4_9mer_WT"),
                       AnnotationFactory.build_annotation(value=self.mhcI_score_allele_9mer_WT, name="bestHLA_allele_netmhcpan4_9mer_Wt"),
                       AnnotationFactory.build_annotation(value=self.mhcI_affinity_9mer_WT, name="best_affinity_netmhcpan4_9mer_WT"),
                       AnnotationFactory.build_annotation(value=self.mhcI_affinity_allele_9mer_WT,
                                  name="bestHLA_allele_affinity_netmhcpan4_9mer_WT"),
                       AnnotationFactory.build_annotation(value=self.mhcI_affinity_epitope_9mer_WT,
                                  name="best_affinity_epitope_netmhcpan4_9mer_WT"),
                       AnnotationFactory.build_annotation(value=self.MHC_epitope_scores, name="MB_epitope_scores"),
                       AnnotationFactory.build_annotation(value=self.MHC_epitope_seqs, name="MB_epitope_sequences"),
                       AnnotationFactory.build_annotation(value=self.MHC_epitope_alleles, name="MB_alleles"),
                       AnnotationFactory.build_annotation(value=self.MHC_number_strong_binders, name="MB_number_pep_MHCscore<1"),
                       AnnotationFactory.build_annotation(value=self.MHC_number_weak_binders, name="MB_number_pep_MHCscore<2"),
                       # generator rate
                       AnnotationFactory.build_annotation(value=self.epitope_affinities, name="MB_affinities"),
                       AnnotationFactory.build_annotation(value=self.generator_rate, name="Generator_rate"),
                       # multiplexed representation WT
                       AnnotationFactory.build_annotation(value=self.MHC_epitope_scores_WT, name="MB_epitope_WT_scores"),
                       AnnotationFactory.build_annotation(value=self.MHC_epitope_seqs_WT, name="MB_epitope_WT_sequences"),
                       AnnotationFactory.build_annotation(value=self.MHC_epitope_alleles_WT, name="MB_alleles_WT"),
                       AnnotationFactory.build_annotation(value=self.MHC_number_strong_binders_WT, name="MB_number_pep_WT_MHCscore<1"),
                       AnnotationFactory.build_annotation(value=self.MHC_number_weak_binders_WT, name="MB_number_pep_WT_MHCscore<2"),
                       AnnotationFactory.build_annotation(value=self.epitope_affinities_WT, name="MB_affinities_WT"),
                       AnnotationFactory.build_annotation(value=self.generator_rate_WT, name="Generator_rate_WT")
                       ]

        for sc, mn in zip(self.MHC_score_all_epitopes, self.mean_type):
            annotations.append(AnnotationFactory.build_annotation(value=sc, name="MB_score_all_epitopes_" + mn))
        for sc, mn in zip(self.MHC_score_top10, self.mean_type):
            annotations.append(AnnotationFactory.build_annotation(value=sc, name="MB_score_top10_" + mn))
        for sc, mn in zip(self.MHC_score_best_per_alelle, self.mean_type):
            # rename MB_score_best_per_alelle_harmonic to PHBR (described in Marty et al)
            annotations.append(
                AnnotationFactory.build_annotation(value=sc, name="MB_score_best_per_alelle_" + mn if mn != 'harmonic' else "PHBR-I"))
        for sc, mn in zip(self.MHC_score_top10_WT, self.mean_type):
            annotations.append(AnnotationFactory.build_annotation(value=sc, name="MB_score_WT_top10_" + mn))
        for sc, mn in zip(self.MHC_score_all_epitopes_WT, self.mean_type):
            annotations.append(AnnotationFactory.build_annotation(value=sc, name="MB_score_WT_all_epitopes_" + mn))
        for sc, mn in zip(self.MHC_score_best_per_alelle_WT, self.mean_type):
            annotations.append(
                AnnotationFactory.build_annotation(value=sc, name="MB_score_WT_best_per_alelle_" + mn if mn != "harmonic" else "PHBR-I_WT"))
        annotations.extend(self._get_positions_and_mutation_in_anchor())
        return annotations

    def _get_positions_and_mutation_in_anchor(self):
        """
        returns if mutation is in anchor position for best affinity epitope over all lengths and best 9mer affinity
        """
        position = EpitopeHelper.position_of_mutation_epitope(
            wild_type=self.best4_affinity_epitope_WT, mutation=self.best4_affinity_epitope)
        position_9mer = EpitopeHelper.position_of_mutation_epitope(
            wild_type=self.mhcI_affinity_epitope_9mer_WT, mutation=self.mhcI_affinity_epitope_9mer)

        return [
            AnnotationFactory.build_annotation(value=position, name="pos_MUT_MHCI_affinity_epi"),
            AnnotationFactory.build_annotation(value=position_9mer, name="pos_MUT_MHCI_affinity_epi_9mer"),
            AnnotationFactory.build_annotation(value=EpitopeHelper.position_of_mutation_epitope(
                wild_type=self.best4_mhc_epitope_WT, mutation=self.best4_mhc_epitope), name="pos_MUT_MHCI_rank_epi"),
            AnnotationFactory.build_annotation(value=EpitopeHelper.position_in_anchor_position(
                position_mhci=position, peptide_length=len(self.best4_affinity_epitope)),
                name="Mutation_in_anchor_netmhcpan"),
            AnnotationFactory.build_annotation(value=EpitopeHelper.position_in_anchor_position(
                position_mhci=position_9mer, peptide_length=9),
                name="Mutation_in_anchor_netmhcpan_9mer"),
            AnnotationFactory.build_annotation(value=EpitopeHelper.position_in_anchor_position(
                position_mhci=position_9mer, peptide_length=len(self.best4_mhc_epitope)),
                name="Mutation_in_anchor_netmhcpan_rank")
            ]
