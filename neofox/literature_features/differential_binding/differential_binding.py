#!/usr/bin/env python

from typing import List

from neofox.model.neoantigen import Annotation
from neofox.model.wrappers import AnnotationFactory
from neofox.predictors.netmhcpan4.combine_netmhcIIpan_pred_multiple_binders import BestAndMultipleBinderMhcII
from neofox.predictors.netmhcpan4.combine_netmhcpan_pred_multiple_binders import BestAndMultipleBinder
from neofox.literature_features.differential_binding.amplitude import Amplitude


class DifferentialBinding:

    def dai(self, score_mutation, score_wild_type, affin_filtering=False):
        """
        Calculates DAI: Returns difference between wt and mut MHC binding score.
        """
        score = None
        try:
            if affin_filtering:
                if score_mutation < 500.0:
                    score = score_wild_type - score_mutation
            else:
                score = score_wild_type - score_mutation
        except TypeError:
            pass
        return score

    def diff_number_binders(self, num_mutation, num_wild_type):
        """
        returns difference of potential candidate epitopes between mutated and wt epitope
        """
        difference = None
        try:
            difference = num_mutation - num_wild_type
        except TypeError:
            pass
        return difference

    def ratio_number_binders(self, num_mutation, num_wild_type):
        """
        returns ratio of number of potential candidate epitopes between mutated and wt epitope.
        if no WT candidate epitopes, returns number of mutated candidate epitopes per mps
        """
        ratio = None
        try:
            ratio = num_mutation / num_wild_type
        except (ZeroDivisionError, TypeError):
            pass
        return ratio

    def classify_adn_cdn(self, score_mutation, amplitude, bdg_cutoff_classical, bdg_cutoff_alternative,
                         amplitude_cutoff,
                         category):
        """
        returns if an epitope belongs to classically and alternatively defined neoepitopes (CDN vs ADN)
        (indicate which category to examine by category)--> Rech et al, 2018
        grouping is based on affinity and affinitiy foldchange between wt and mut
        """
        group = None
        try:
            if category == "CDN":
                group = score_mutation < bdg_cutoff_classical
            elif category == "ADN":
                group = score_mutation < bdg_cutoff_alternative and amplitude > amplitude_cutoff
        except (ValueError, TypeError):
            pass
        return group

    def get_annotations_dai(self, netmhcpan: BestAndMultipleBinder, netmhc2pan: BestAndMultipleBinderMhcII) -> List[Annotation]:
        return [
            # MHC I
            AnnotationFactory.build_annotation(
                name="DAI_MHCI_affinity_cutoff500nM", value=self.dai(
                    score_mutation=netmhcpan.best4_affinity, score_wild_type=netmhcpan.best4_affinity_WT,
                    affin_filtering=True)),
            AnnotationFactory.build_annotation(
                name="DAI_MHCI_affinity", value=self.dai(
                    score_mutation=netmhcpan.best4_affinity, score_wild_type=netmhcpan.best4_affinity_WT)),
            AnnotationFactory.build_annotation(
                name="DAI_MHCI_rank", value=self.dai(
                    score_mutation=netmhcpan.best4_mhc_score, score_wild_type=netmhcpan.best4_mhc_score_WT)),
            AnnotationFactory.build_annotation(
                name="DAI_MHCI_multiple_binding", value=self.dai(
                    score_mutation=netmhcpan.MHC_score_top10[1],
                    score_wild_type=netmhcpan.MHC_score_top10_WT[1])),
            # MHC II
            AnnotationFactory.build_annotation(
                value=self.dai(score_mutation=netmhc2pan.best_mhcII_pan_affinity,
                               score_wild_type=netmhc2pan.best_mhcII_affinity_WT, affin_filtering=True),
                name="DAI_MHCII_affinity_cutoff500nM"),
            AnnotationFactory.build_annotation(
                value=self.dai(score_mutation=netmhc2pan.best_mhcII_pan_affinity,
                               score_wild_type=netmhc2pan.best_mhcII_affinity_WT),
                name="DAI_MHCII_affinity"),
            AnnotationFactory.build_annotation(
                value=self.dai(score_mutation=netmhc2pan.best_mhcII_pan_score,
                               score_wild_type=netmhc2pan.best_mhcII_pan_score_WT),
                name="DAI_MHCII_rank")
        ]


    def get_annotations(self, netmhcpan: BestAndMultipleBinder, amplitude: Amplitude) -> List[Annotation]:

        bdg_cutoff_classical_mhci = 50
        bdg_cutoff_alternative_mhci = 5000
        amplitude_cutoff_mhci = 10

        return [
            AnnotationFactory.build_annotation(name="CDN_MHCI", value=self.classify_adn_cdn(
                score_mutation=netmhcpan.best4_affinity, amplitude=amplitude.amplitude_mhci_affinity,
                bdg_cutoff_classical=bdg_cutoff_classical_mhci, bdg_cutoff_alternative=bdg_cutoff_alternative_mhci,
                amplitude_cutoff=amplitude_cutoff_mhci, category="CDN")),
            AnnotationFactory.build_annotation(name="ADN_MHCI", value=self.classify_adn_cdn(
                score_mutation=netmhcpan.best4_affinity, amplitude=amplitude.amplitude_mhci_affinity,
                bdg_cutoff_classical=bdg_cutoff_classical_mhci, bdg_cutoff_alternative=bdg_cutoff_alternative_mhci,
                amplitude_cutoff=amplitude_cutoff_mhci, category="ADN")),
            AnnotationFactory.build_annotation(value=self.diff_number_binders(
                num_mutation=netmhcpan.MHC_number_strong_binders, num_wild_type=netmhcpan.MHC_number_strong_binders_WT),
                name="Difference_number_epitopes_MHCI_strong_binder"),
            AnnotationFactory.build_annotation(value=self.diff_number_binders(
                num_mutation=netmhcpan.MHC_number_weak_binders, num_wild_type=netmhcpan.MHC_number_strong_binders_WT),
                name="Difference_number_epitopes_MHCI_weak_binder"),
            AnnotationFactory.build_annotation(value=self.ratio_number_binders(
                num_mutation=netmhcpan.MHC_number_strong_binders, num_wild_type=netmhcpan.MHC_number_strong_binders_WT),
                name="Ratio_number_epitopes_MHCI_strong_binder}"),
            AnnotationFactory.build_annotation(value=self.ratio_number_binders(
                num_mutation=netmhcpan.MHC_number_weak_binders, num_wild_type=netmhcpan.MHC_number_weak_binders_WT),
                name="Ratio_number_epitopes_MHCI_weak_binder")
        ]

    def get_annotations_mhc2(self, netmhc2pan: BestAndMultipleBinderMhcII, amplitude: Amplitude) -> List[Annotation]:

        bdg_cutoff_classical_mhcii = 1
        bdg_cutoff_alternative_mhcii = 4
        amplitude_cutoff_mhcii = 4

        return [
            AnnotationFactory.build_annotation(
                value=self.classify_adn_cdn(
                    score_mutation=netmhc2pan.best_mhcII_pan_score, amplitude=amplitude.amplitude_mhcii_rank,
                    bdg_cutoff_classical=bdg_cutoff_classical_mhcii,
                    bdg_cutoff_alternative=bdg_cutoff_alternative_mhcii, amplitude_cutoff=amplitude_cutoff_mhcii,
                    category="CDN"),
                name="CDN_MHCII"),
            AnnotationFactory.build_annotation(
                value=self.classify_adn_cdn(
                    score_mutation=netmhc2pan.best_mhcII_pan_score, amplitude=amplitude.amplitude_mhcii_rank,
                    bdg_cutoff_classical=bdg_cutoff_classical_mhcii,
                    bdg_cutoff_alternative=bdg_cutoff_alternative_mhcii, amplitude_cutoff=amplitude_cutoff_mhcii,
                    category="ADN"),
                name="ADN_MHCII"),
            AnnotationFactory.build_annotation(value=self.diff_number_binders(
                num_mutation=netmhc2pan.MHCII_number_strong_binders,
                num_wild_type=netmhc2pan.MHCII_number_strong_binders_WT),
                name="Difference_number_epitopes_MHCII_strong_binder"),
            AnnotationFactory.build_annotation(value=self.diff_number_binders(
                num_mutation=netmhc2pan.MHCII_number_weak_binders,
                num_wild_type=netmhc2pan.MHCII_number_weak_binders_WT),
                name="Difference_number_epitopes_MHCII_weak_binder"),
            AnnotationFactory.build_annotation(value=self.ratio_number_binders(
                num_mutation=netmhc2pan.MHCII_number_strong_binders,
                num_wild_type=netmhc2pan.MHCII_number_strong_binders_WT),
                name="Ratio_number_epitopes_MHCII_strong_binder"),
            AnnotationFactory.build_annotation(value=self.ratio_number_binders(
                num_mutation=netmhc2pan.MHCII_number_weak_binders,
                num_wild_type=netmhc2pan.MHCII_number_weak_binders_WT),
                name="Ratio_number_epitopes_MHCII_weak_binder")
        ]
