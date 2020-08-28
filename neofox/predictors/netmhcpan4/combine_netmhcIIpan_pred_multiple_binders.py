#!/usr/bin/env python
from typing import List

import neofox.predictors.netmhcpan4.netmhcIIpan_prediction as netmhcIIpan_prediction
from neofox import MHC_II
from neofox.helpers import intermediate_files
from neofox.literature_features.differential_binding import DifferentialBinding
from neofox.model.neoantigen import Annotation
from neofox.model.wrappers import AnnotationFactory
from neofox.predictors.netmhcpan4 import multiple_binders


class BestAndMultipleBinderMhcII:

    def __init__(self, runner, configuration):
        """
        :type runner: neofox.helpers.runner.Runner
        :type configuration: neofox.references.DependenciesConfiguration
        """
        self.runner = runner
        self.configuration = configuration
        self.differential_binding = DifferentialBinding()
        self.mean_type = ["arithmetic", "harmonic", "geometric"]
        self.MHCII_score_all_epitopes = ["NA", "NA", "NA"]
        self.MHCII_score_top10 = ["NA", "NA", "NA"]
        self.MHCII_score_best_per_alelle = ["NA", "NA", "NA"]
        self.MHCII_number_strong_binders = "NA"
        self.MHCII_number_weak_binders = "NA"
        self.MHCII_epitope_seqs = "NA"
        self.MHCII_epitope_scores = "NA"
        self.MHCII_epitope_alleles = "NA"
        self.best_mhcII_pan_score = "NA"
        self.best_mhcII_pan_epitope = "-"
        self.best_mhcII_pan_allele = "NA"
        self.best_mhcII_pan_position = "NA"
        self.best_mhcII_pan_affinity = "NA"
        self.best_mhcII_pan_affinity_epitope = "-"
        self.best_mhcII_pan_affinity_allele = "NA"
        self.best_mhcII_pan_affinity_position = "NA"
        # WT features
        self.MHCII_epitope_scores_WT = "Na"
        self.MHCII_epitope_seqs_WT = "NA"
        self.MHCII_epitope_alleles_WT = "NA"
        self.MHCII_score_top10_WT = ["NA", "NA", "NA"]
        self.MHCII_score_all_epitopes_WT = ["NA", "NA", "NA"]
        self.MHCII_score_best_per_alelle_WT = ["NA", "NA", "NA"]
        self.MHCII_number_strong_binders_WT = "NA"
        self.MHCII_number_weak_binders_WT = "NA"
        self.best_mhcII_pan_score_WT = "NA"
        self.best_mhcII_pan_epitope_WT = "-"
        self.best_mhcII_pan_allele_WT = "NA"
        self.best_mhcII_affinity_WT = "NA"
        self.best_mhcII_affinity_epitope_WT = "-"
        self.best_mhcII_affinity_allele_WT = "NA"

    def MHCII_MB_score_best_per_allele(self, tuple_best_per_allele):
        '''returns list of multiple binding scores for mhcII considering best epitope per allele, applying different types of means (harmonic ==> PHRB-II, Marty et al).
        2 copies of DRA - DRB1 --> consider this gene 2x when averaging mhcii binding scores
        '''
        number_alleles = len(tuple_best_per_allele)
        multbind = multiple_binders.MultipleBinding()
        tuple_best_per_allele_new = list(tuple_best_per_allele)
        for best_epi in tuple_best_per_allele:
            if best_epi[-1].startswith("DRB1"):
                tuple_best_per_allele_new.append(best_epi)
        if len(tuple_best_per_allele_new) == 12:
            # 12 genes gene copies should be included into PHBR_II
            best_scores_allele = multbind.scores_to_list(tuple_best_per_allele_new)
            return multbind.get_means(best_scores_allele)
        else:
            return ["NA", "NA", "NA"]

    def run(self, sequence, sequence_reference, alleles, set_available_mhc):
        '''predicts MHC epitopes; returns on one hand best binder and on the other hand multiple binder analysis is performed
        '''
        ### PREDICTION FOR MUTATED SEQUENCE
        tmp_prediction = intermediate_files.create_temp_file(prefix="netmhcpanpred_", suffix=".csv")
        np = netmhcIIpan_prediction.NetMhcIIPanPredictor(runner=self.runner, configuration=self.configuration)
        mb = multiple_binders.MultipleBinding()
        tmp_fasta = intermediate_files.create_temp_fasta([sequence], prefix="tmp_singleseq_")
        alleles_formated = np.generate_mhcII_alelles_combination_list(alleles, set_available_mhc)
        np.mhcII_prediction(alleles, set_available_mhc, tmp_fasta, tmp_prediction)
        position_xmer_sequence = np.mut_position_xmer_seq(xmer_wt=sequence_reference, xmer_mut=sequence)
        try:
            preds = np.filter_binding_predictions(position_xmer_sequence, tmp_prediction)
            # multiple binding
            list_tups = mb.generate_epi_tuple(preds, mhc=MHC_II)
            self.MHCII_epitope_scores = "/".join([tup[0] for tup in list_tups])
            self.MHCII_epitope_seqs = "/".join([tup[2] for tup in list_tups])
            self.MHCII_epitope_alleles = "/".join([tup[3] for tup in list_tups])
            top10 = mb.extract_top10_epis(list_tups)
            best_per_alelle = mb.extract_best_epi_per_alelle(list_tups, alleles_formated)
            all = mb.scores_to_list(list_tups)
            top10 = mb.scores_to_list(top10)
            self.MHCII_score_top10 = mb.get_means(top10)
            self.MHCII_score_all_epitopes = mb.get_means(all)
            self.MHCII_score_best_per_alelle = self.MHCII_MB_score_best_per_allele(best_per_alelle)
            self.MHCII_number_strong_binders = mb.determine_number_of_binders(all, 2)
            self.MHCII_number_weak_binders = mb.determine_number_of_binders(all, 10)
            # best prediction
            best_epi = np.minimal_binding_score(preds)
            self.best_mhcII_pan_score = np.add_best_epitope_info(best_epi, "%Rank")
            self.best_mhcII_pan_epitope = np.add_best_epitope_info(best_epi, "Peptide")
            self.best_mhcII_pan_allele = np.add_best_epitope_info(best_epi, "Allele")
            self.best_mhcII_pan_position = np.add_best_epitope_info(best_epi, "Seq")
            best_epi_affinity = np.minimal_binding_score(preds, rank=False)
            self.best_mhcII_pan_affinity = np.add_best_epitope_info(best_epi_affinity, "Affinity(nM)")
            self.best_mhcII_pan_affinity_epitope = np.add_best_epitope_info(best_epi_affinity, "Peptide")
            self.best_mhcII_pan_affinity_allele = np.add_best_epitope_info(best_epi_affinity, "Allele")
            self.best_mhcII_pan_affinity_position = np.add_best_epitope_info(best_epi_affinity, "Seq")
        except IndexError:
            # if neofox sequence shorter than 15 aa
            pass

        ### PREDICTION FOR WT SEQUENCE
        tmp_prediction = intermediate_files.create_temp_file(prefix="netmhcpanpred_", suffix=".csv")
        np = netmhcIIpan_prediction.NetMhcIIPanPredictor(runner=self.runner, configuration=self.configuration)
        mb = multiple_binders.MultipleBinding()
        tmp_fasta = intermediate_files.create_temp_fasta([sequence_reference], prefix="tmp_singleseq_")
        np.mhcII_prediction(alleles, set_available_mhc, tmp_fasta, tmp_prediction)
        try:
            preds = np.filter_binding_predictions(position_xmer_sequence, tmp_prediction)
            # multiple binding
            list_tups = mb.generate_epi_tuple(preds, mhc=MHC_II)
            self.MHCII_epitope_scores_WT = "/".join([tup[0] for tup in list_tups])
            self.epitope_affinities__mhcII_pan_WT = "/".join([tup[1] for tup in list_tups])
            self.MHCII_epitope_seqs_WT = "/".join([tup[2] for tup in list_tups])
            self.MHCII_epitope_alleles_WT = "/".join([tup[3] for tup in list_tups])
            top10 = mb.extract_top10_epis(list_tups)
            best_per_alelle = mb.extract_best_epi_per_alelle(list_tups, alleles_formated)
            all = mb.scores_to_list(list_tups)
            all_affinities = mb.affinities_to_list(list_tups)
            top10 = mb.scores_to_list(top10)
            self.MHCII_score_top10_WT = mb.get_means(top10)
            self.MHCII_score_all_epitopes_WT = mb.get_means(all)
            self.MHCII_score_best_per_alelle_WT = self.MHCII_MB_score_best_per_allele(best_per_alelle)
            self.MHCII_number_strong_binders_WT = mb.determine_number_of_binders(all, 1)
            self.MHCII_number_weak_binders_WT = mb.determine_number_of_binders(all, 2)
            # best prediction
            best_epi = np.filter_for_wt_epitope_position(preds, self.best_mhcII_pan_epitope,
                                                         position_epi_xmer=self.best_mhcII_pan_position)
            self.best_mhcII_pan_score_WT = np.add_best_epitope_info(best_epi, "%Rank")
            self.best_mhcII_pan_epitope_WT = np.add_best_epitope_info(best_epi, "Peptide")
            self.best_mhcII_pan_allele_WT = np.add_best_epitope_info(best_epi, "Allele")
            best_epi_affinity = np.filter_for_wt_epitope_position(preds, self.best_mhcII_pan_affinity_epitope,
                                                                  position_epi_xmer=self.best_mhcII_pan_affinity_position)
            self.best_mhcII_affinity_WT = np.add_best_epitope_info(best_epi_affinity, "Affinity(nM)")
            self.best_mhcII_affinity_epitope_WT = np.add_best_epitope_info(best_epi_affinity, "Peptide")
            self.best_mhcII_affinity_allele_WT = np.add_best_epitope_info(best_epi_affinity, "Allele")
        except IndexError:
            # if neofox sequence shorter than 15 aa
            pass

    def get_annotations(self) -> List[Annotation]:
        annotations =  [
            # netmhcpan4 MUT scores
            AnnotationFactory.build_annotation(value=self.best_mhcII_pan_score, name="best%Rank_netmhcIIpan"),
            AnnotationFactory.build_annotation(value=self.best_mhcII_pan_epitope, name="best_epitope_netmhcIIpan"),
            AnnotationFactory.build_annotation(value=self.best_mhcII_pan_allele, name="bestHLA_allele_netmhcIIpan"),
            # netmhcpan4 mut affinity
            AnnotationFactory.build_annotation(value=self.best_mhcII_pan_affinity, name="best_affinity_netmhcIIpan"),
            AnnotationFactory.build_annotation(value=self.best_mhcII_pan_affinity_epitope,
                                               name="best_affinity_epitope_netmhcIIpan"),
            AnnotationFactory.build_annotation(value=self.best_mhcII_pan_affinity_allele,
                                               name="bestHLA_allele_affinity_netmhcIIpan"),
            # netmhcIIpan WT scores
            AnnotationFactory.build_annotation(value=self.best_mhcII_pan_score_WT, name="best%Rank_netmhcIIpan_WT"),
            AnnotationFactory.build_annotation(value=self.best_mhcII_pan_epitope_WT, name="best_epitope_netmhcIIpan_WT"),
            AnnotationFactory.build_annotation(value=self.best_mhcII_pan_allele_WT, name="bestHLA_allele_netmhcIIpan_Wt"),
            # netmhcIIpan wt affinity
            AnnotationFactory.build_annotation(value=self.best_mhcII_affinity_WT, name="best_affinity_netmhcIIpan_WT"),
            AnnotationFactory.build_annotation(value=self.best_mhcII_affinity_epitope_WT,
                                               name="best_affinity_epitope_netmhcIIpan_WT"),
            AnnotationFactory.build_annotation(value=self.best_mhcII_affinity_allele_WT,
                                               name="bestHLA_allele_affinity_netmhcIIpan_WT"),

            AnnotationFactory.build_annotation(value=self.MHCII_epitope_scores, name="MB_mhcII_epitope_scores"),
            AnnotationFactory.build_annotation(value=self.MHCII_epitope_seqs, name="MB_mhcII_epitope_sequences"),
            AnnotationFactory.build_annotation(value=self.MHCII_epitope_alleles, name="MB_mhcII_alleles"),
            AnnotationFactory.build_annotation(value=self.MHCII_number_strong_binders,
                                               name="MB_number_pep_MHCIIscore<2"),
            AnnotationFactory.build_annotation(value=self.MHCII_number_weak_binders,
                                               name="MB_number_pep_MHCIIscore<10"),

            AnnotationFactory.build_annotation(value=self.MHCII_epitope_scores_WT, name="MB_mhcII_epitope_scores_WT"),
            AnnotationFactory.build_annotation(value=self.MHCII_epitope_seqs_WT, name="MB_mhcII_epitope_sequences_WT"),
            AnnotationFactory.build_annotation(value=self.MHCII_epitope_alleles_WT, name="MB_mhcII_alleles_WT"),
            AnnotationFactory.build_annotation(value=self.MHCII_number_strong_binders_WT,
                                               name="MB_number_pep_MHCIIscore<2_WT"),
            AnnotationFactory.build_annotation(value=self.MHCII_number_weak_binders_WT,
                                               name="MB_number_pep_MHCIIscore<10_WT")
            ]
        # multiplexed representation MUT MHC II
        for sc, mn in zip(self.MHCII_score_all_epitopes, self.mean_type):
            annotations.append(AnnotationFactory.build_annotation(value=sc, name="MB_score_MHCII_all_epitopes_" + mn))
        for sc, mn in zip(self.MHCII_score_top10, self.mean_type):
            annotations.append(AnnotationFactory.build_annotation(value=sc, name="MB_score_MHCII_top10_" + mn))
        for sc, mn in zip(self.MHCII_score_best_per_alelle, self.mean_type):
            # rename MB_score_best_per_alelle_harmonic to PHBR (described in Marty et al)
            annotations.append(AnnotationFactory.build_annotation(
                value=sc, name="MB_score_MHCII_best_per_alelle_" + mn if mn != "harmonic" else "PHBR-II"))
        # multiplexed representation WT MHC II
        for sc, mn in zip(self.MHCII_score_all_epitopes_WT, self.mean_type):
            annotations.append(AnnotationFactory.build_annotation(value=sc,
                                                                  name="MB_score_MHCII_all_epitopes_WT_" + mn))
        for sc, mn in zip(self.MHCII_score_top10_WT, self.mean_type):
            annotations.append(AnnotationFactory.build_annotation(value=sc, name="MB_score_MHCII_top10_WT_" + mn))
        for sc, mn in zip(self.MHCII_score_best_per_alelle_WT, self.mean_type):
            # rename MB_score_best_per_alelle_harmonic to PHBR (described in Marty et al)
            annotations.append(AnnotationFactory.build_annotation(
                value=sc, name="MB_score_MHCII_best_per_alelle_WT_" + mn if mn != "harmonic" else "PHBR-II_WT"))

        annotations.extend(self.get_differential_binding_annotations(
            aff_mut=self.best_mhcII_pan_affinity, aff_wt=self.best_mhcII_affinity_WT,
            rank_mut=self.best_mhcII_pan_score, rank_wt=self.best_mhcII_pan_score_WT
        ))
        annotations.extend(self.get_multiple_binding_annotations())

        return annotations

    def get_differential_binding_annotations(self, aff_mut, aff_wt, rank_mut, rank_wt) -> List[Annotation]:
        """
        returns DAI for MHC II based on affinity (filtered + no filtered) and rank
        """
        # dai mhc II affinity
        return [
            AnnotationFactory.build_annotation(
                value=self.differential_binding.dai(score_mutation=aff_mut, score_wild_type=aff_wt),
                name="DAI_mhcII_affinity"),
            AnnotationFactory.build_annotation(
                value=self.differential_binding.dai(score_mutation=aff_mut, score_wild_type=aff_wt, affin_filtering=True),
                name="DAI_mhcII_affinity_aff_filtered"),
            # dai mhc II netMHCIIpan score
            AnnotationFactory.build_annotation(
                value=self.differential_binding.dai(score_mutation=rank_mut, score_wild_type=rank_wt),
                name="DAI_mhcII_rank")
        ]

    def get_multiple_binding_annotations(self):

        num_strong_binders_mutation = self.MHCII_number_strong_binders
        num_strong_binders_wild_type = self.MHCII_number_strong_binders_WT
        num_weak_binders_mutation = self.MHCII_number_weak_binders
        num_weak_binders_wild_type = self.MHCII_number_weak_binders_WT
        return [
            AnnotationFactory.build_annotation(value=self.differential_binding.diff_number_binders(
                num_mutation=num_strong_binders_mutation, num_wild_type=num_strong_binders_wild_type),
                name="Diff_numb_epis_mhcII<2"),
            AnnotationFactory.build_annotation(value=self.differential_binding.ratio_number_binders(
                num_mutation=num_strong_binders_mutation, num_wild_type=num_strong_binders_wild_type),
                name="Ratio_numb_epis_mhcII<2"),
            AnnotationFactory.build_annotation(value=self.differential_binding.diff_number_binders(
                num_mutation=num_weak_binders_mutation, num_wild_type=num_weak_binders_wild_type),
                name="Diff_numb_epis_mhcII<10"),
            AnnotationFactory.build_annotation(value=self.differential_binding.ratio_number_binders(
                num_mutation=num_weak_binders_mutation, num_wild_type=num_weak_binders_wild_type),
                name="Ratio_numb_epis_mhcII<10")
        ]
