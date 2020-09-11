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
from typing import List

import neofox.MHC_predictors.netmhcpan.netmhcIIpan_prediction as netmhcIIpan_prediction
from neofox import MHC_II
from neofox.helpers import intermediate_files
from neofox.model.neoantigen import Annotation
from neofox.model.wrappers import AnnotationFactory
from neofox.MHC_predictors.netmhcpan import multiple_binders
import neofox.helpers.casting as casting


class BestAndMultipleBinderMhcII:

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
        self.MHCII_score_all_epitopes = [None, None, None]
        self.MHCII_score_top10 = [None, None, None]
        self.MHCII_score_best_per_alelle = [None, None, None]
        self.MHCII_number_strong_binders = None
        self.MHCII_number_weak_binders = None
        self.MHCII_epitope_seqs = None
        self.MHCII_epitope_scores = None
        self.MHCII_epitope_alleles = None
        self.best_mhcII_pan_score = None
        self.best_mhcII_pan_epitope = "-"
        self.best_mhcII_pan_allele = None
        self.best_mhcII_pan_position = None
        self.best_mhcII_pan_affinity = None
        self.best_mhcII_pan_affinity_epitope = "-"
        self.best_mhcII_pan_affinity_allele = None
        self.best_mhcII_pan_affinity_position = None
        # WT features
        self.MHCII_epitope_scores_WT = None
        self.MHCII_epitope_seqs_WT = None
        self.MHCII_epitope_alleles_WT = None
        self.MHCII_score_top10_WT = [None, None, None]
        self.MHCII_score_all_epitopes_WT = [None, None, None]
        self.MHCII_score_best_per_alelle_WT = [None, None, None]
        self.MHCII_number_strong_binders_WT = None
        self.MHCII_number_weak_binders_WT = None
        self.best_mhcII_pan_score_WT = None
        self.best_mhcII_pan_epitope_WT = "-"
        self.best_mhcII_pan_allele_WT = None
        self.best_mhcII_affinity_WT = None
        self.best_mhcII_affinity_epitope_WT = "-"
        self.best_mhcII_affinity_allele_WT = None

    def MHCII_MB_score_best_per_allele(self, tuple_best_per_allele):
        """returns list of multiple binding scores for mhcII considering best epitope per allele,
        applying different types of means (harmonic ==> PHRB-II, Marty et al).
        2 copies of DRA - DRB1 --> consider this gene 2x when averaging mhcii binding scores
        """
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
            return [None, None, None]

    def run(self, sequence, sequence_reference, alleles, set_available_mhc):
        """predicts MHC epitopes; returns on one hand best binder and on the other hand multiple binder analysis is performed
        """
        # mutation
        self._initialise()
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
            self.MHCII_epitope_scores = "/".join([str(tup[0]) for tup in list_tups])
            self.MHCII_epitope_seqs = "/".join([str(tup[2]) for tup in list_tups])
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
            self.best_mhcII_pan_score = casting.to_float(np.add_best_epitope_info(best_epi, "%Rank"))
            self.best_mhcII_pan_epitope = np.add_best_epitope_info(best_epi, "Peptide")
            self.best_mhcII_pan_allele = np.add_best_epitope_info(best_epi, "Allele")
            self.best_mhcII_pan_position = np.add_best_epitope_info(best_epi, "Seq")
            best_epi_affinity = np.minimal_binding_score(preds, rank=False)
            self.best_mhcII_pan_affinity = casting.to_float(np.add_best_epitope_info(best_epi_affinity, "Affinity(nM)"))
            self.best_mhcII_pan_affinity_epitope = np.add_best_epitope_info(best_epi_affinity, "Peptide")
            self.best_mhcII_pan_affinity_allele = np.add_best_epitope_info(best_epi_affinity, "Allele")
            self.best_mhcII_pan_affinity_position = np.add_best_epitope_info(best_epi_affinity, "Seq")
        except IndexError:
            # if neofox sequence shorter than 15 aa
            pass

        # wt
        tmp_prediction = intermediate_files.create_temp_file(prefix="netmhcpanpred_", suffix=".csv")
        np = netmhcIIpan_prediction.NetMhcIIPanPredictor(runner=self.runner, configuration=self.configuration)
        mb = multiple_binders.MultipleBinding()
        tmp_fasta = intermediate_files.create_temp_fasta([sequence_reference], prefix="tmp_singleseq_")
        np.mhcII_prediction(alleles, set_available_mhc, tmp_fasta, tmp_prediction)
        try:
            preds = np.filter_binding_predictions(position_xmer_sequence, tmp_prediction)
            # multiple binding
            list_tups = mb.generate_epi_tuple(preds, mhc=MHC_II)
            self.MHCII_epitope_scores_WT = "/".join([str(tup[0]) for tup in list_tups])
            self.epitope_affinities__mhcII_pan_WT = "/".join([str(tup[1]) for tup in list_tups])
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
            self.best_mhcII_pan_score_WT = casting.to_float(np.add_best_epitope_info(best_epi, "%Rank"))
            self.best_mhcII_pan_epitope_WT = np.add_best_epitope_info(best_epi, "Peptide")
            self.best_mhcII_pan_allele_WT = np.add_best_epitope_info(best_epi, "Allele")
            best_epi_affinity = np.filter_for_wt_epitope_position(preds, self.best_mhcII_pan_affinity_epitope,
                                                                  position_epi_xmer=self.best_mhcII_pan_affinity_position)
            self.best_mhcII_affinity_WT = casting.to_float(np.add_best_epitope_info(best_epi_affinity, "Affinity(nM)"))
            self.best_mhcII_affinity_epitope_WT = np.add_best_epitope_info(best_epi_affinity, "Peptide")
            self.best_mhcII_affinity_allele_WT = np.add_best_epitope_info(best_epi_affinity, "Allele")
        except IndexError:
            # if neofox sequence shorter than 15 aa
            pass

    def get_annotations(self) -> List[Annotation]:
        annotations =  [
            AnnotationFactory.build_annotation(value=self.best_mhcII_pan_score, name="Best_rank_MHCII_score"),
            AnnotationFactory.build_annotation(value=self.best_mhcII_pan_epitope, name="Best_rank_MHCII_score_epitope"),
            AnnotationFactory.build_annotation(value=self.best_mhcII_pan_allele, name="Best_rank_MHCII_score_allele"),
            AnnotationFactory.build_annotation(value=self.best_mhcII_pan_affinity, name="Best_affinity_MHCII_score"),
            AnnotationFactory.build_annotation(value=self.best_mhcII_pan_affinity_epitope,
                                               name="Best_affinity_MHCII_epitope"),
            AnnotationFactory.build_annotation(value=self.best_mhcII_pan_affinity_allele,
                                               name="Best_affinity_MHCII_allele"),
            AnnotationFactory.build_annotation(value=self.best_mhcII_pan_score_WT, name="Best_rank_MHCII_score_WT"),
            AnnotationFactory.build_annotation(value=self.best_mhcII_pan_epitope_WT,
                                               name="Best_rank_MHCII_score_epitope_WT"),
            AnnotationFactory.build_annotation(value=self.best_mhcII_pan_allele_WT,
                                               name="Best_rank_MHCII_score_allele_WT"),
            AnnotationFactory.build_annotation(value=self.best_mhcII_affinity_WT, name="Best_affinity_MHCII_score_WT"),
            AnnotationFactory.build_annotation(value=self.best_mhcII_affinity_epitope_WT,
                                               name="Best_affinity_MHCII_epitope_WT"),
            AnnotationFactory.build_annotation(value=self.best_mhcII_affinity_allele_WT,
                                               name="Best_affinity_MHCII_allele_WT"),

            AnnotationFactory.build_annotation(value=self.MHCII_epitope_scores, name="All_ranks_MHCII"),
            AnnotationFactory.build_annotation(value=self.MHCII_epitope_seqs, name="All_epitopes_MHCII"),
            AnnotationFactory.build_annotation(value=self.MHCII_epitope_alleles, name="All_alleles_MHCII"),
            AnnotationFactory.build_annotation(value=self.MHCII_number_strong_binders,
                                               name="Number_strong_binders_MHCII"),
            AnnotationFactory.build_annotation(value=self.MHCII_number_weak_binders,
                                               name="Number_weak_binders_MHCII"),

            AnnotationFactory.build_annotation(value=self.MHCII_epitope_scores_WT, name="All_ranks_MHCII_WT"),
            AnnotationFactory.build_annotation(value=self.MHCII_epitope_seqs_WT, name="All_epitopes_MHCII_WT"),
            AnnotationFactory.build_annotation(value=self.MHCII_epitope_alleles_WT, name="All_alleles_MHCII_WT"),
            AnnotationFactory.build_annotation(value=self.MHCII_number_strong_binders_WT,
                                               name="Number_strong_binders_MHCII_WT"),
            AnnotationFactory.build_annotation(value=self.MHCII_number_weak_binders_WT,
                                               name="Number_weak_binders_MHCII_WT")
            ]
        # multiplexed representation MUT MHC II
        for sc, mn in zip(self.MHCII_score_all_epitopes, self.mean_type):
            annotations.append(AnnotationFactory.build_annotation(value=sc,
                                                                  name="Multiple_binding_score_MHCII_all_epitopes_" + mn))
        for sc, mn in zip(self.MHCII_score_top10, self.mean_type):
            annotations.append(AnnotationFactory.build_annotation(value=sc, name="Multiple_binding_score_MHCII_top10_" + mn))
        for sc, mn in zip(self.MHCII_score_best_per_alelle, self.mean_type):
            # rename MB_score_best_per_alelle_harmonic to PHBR (described in Marty et al)
            annotations.append(AnnotationFactory.build_annotation(
                value=sc, name="Multiple_binding_score_MHCII_best_per_alelle_" + mn if mn != "harmonic" else "PHBR-II"))
        # multiplexed representation WT MHC II
        for sc, mn in zip(self.MHCII_score_all_epitopes_WT, self.mean_type):
            annotations.append(AnnotationFactory.build_annotation(value=sc,
                                                                  name="Multiple_binding_score_MHCII_all_epitopes_WT_" + mn))
        for sc, mn in zip(self.MHCII_score_top10_WT, self.mean_type):
            annotations.append(AnnotationFactory.build_annotation(value=sc,
                                                                  name="Multiple_binding_score_MHCII_top10_WT_" + mn))
        for sc, mn in zip(self.MHCII_score_best_per_alelle_WT, self.mean_type):
            # rename MB_score_best_per_alelle_harmonic to PHBR (described in Marty et al)
            annotations.append(AnnotationFactory.build_annotation(
                value=sc, name="Multiple_binding_score_MHCII_best_per_alelle_WT_" + mn if mn != "harmonic" else "PHBR-II_WT"))

        return annotations
