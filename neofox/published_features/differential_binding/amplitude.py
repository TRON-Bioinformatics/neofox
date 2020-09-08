#!/usr/bin/env python

from typing import List

from neofox.model.neoantigen import Annotation
from neofox.model.wrappers import AnnotationFactory
from neofox.MHC_predictors.netmhcpan.combine_netmhcIIpan_pred_multiple_binders import BestAndMultipleBinderMhcII
from neofox.MHC_predictors.netmhcpan.combine_netmhcpan_pred_multiple_binders import BestAndMultipleBinder


class Amplitude:

    def __init__(self):
        self.amplitude_mhci_affinity = None
        self.amplitude_mhci_affinity_9mer = None
        self.amplitude_mhci_rank = None
        self.amplitude_mhci_MB = None
        self.amplitude_mhcii_affinity = None
        self.amplitude_mhcii_rank = None
        self.amplitude_mhcii_MB = None

    def calculate_amplitude_mhc(self, score_mutation, score_wild_type, apply_correction=False):
        """
        This function calculates the amplitude between mutated and wt epitope according to Balachandran et al.
        when affinity is used, use correction from Luksza et al. *1/(1+0.0003*aff_wt)
        """
        amplitude_mhc = None
        try:
            candidate_amplitude_mhc = score_wild_type / score_mutation
            if apply_correction:  # nine_mer or affinity:
                amplitude_mhc = candidate_amplitude_mhc * self._calculate_correction(score_wild_type)
            else:
                amplitude_mhc = candidate_amplitude_mhc
        except(ZeroDivisionError, ValueError, TypeError):
            pass
        return amplitude_mhc

    def _calculate_correction(self, score_wild_type):
        return 1 / (1 + 0.0003 * float(score_wild_type))

    def run(self, netmhcpan: BestAndMultipleBinder, netmhc2pan: BestAndMultipleBinderMhcII):
        # MHC I
        self.amplitude_mhci_affinity = self.calculate_amplitude_mhc(
            score_mutation=netmhcpan.best4_affinity, score_wild_type=netmhcpan.best4_affinity_WT,
            apply_correction=True)
        self.amplitude_mhci_rank = self.calculate_amplitude_mhc(
            score_mutation=netmhcpan.best4_mhc_score, score_wild_type=netmhcpan.best4_mhc_score_WT)
        self.amplitude_mhci_affinity_9mer = self.calculate_amplitude_mhc(
            score_mutation=netmhcpan.mhcI_affinity_9mer, score_wild_type=netmhcpan.mhcI_affinity_9mer_WT,
            apply_correction=True)
        self.amplitude_mhci_MB = self.calculate_amplitude_mhc(
            score_mutation=netmhcpan.MHC_score_top10[1],
            score_wild_type=netmhcpan.MHC_score_top10_WT[1])
        # MHC II
        self.amplitude_mhcii_affinity = self.calculate_amplitude_mhc(
            score_mutation=netmhc2pan.best_mhcII_pan_affinity, score_wild_type=netmhc2pan.best_mhcII_affinity_WT,
            apply_correction=True)
        self.amplitude_mhcii_rank = self.calculate_amplitude_mhc(
            score_mutation=netmhc2pan.best_mhcII_pan_score, score_wild_type=netmhc2pan.best_mhcII_pan_score_WT)
        self.amplitude_mhcii_MB = self.calculate_amplitude_mhc(
            score_mutation=netmhc2pan.MHCII_score_top10[1], score_wild_type=netmhc2pan.MHCII_score_top10_WT[1])

    def get_annotations(self) -> List[Annotation]:
        return [
            AnnotationFactory.build_annotation(
                value=self.amplitude_mhci_affinity, name="Amplitude_MHCI_affinity"),
            AnnotationFactory.build_annotation(value=self.amplitude_mhci_rank, name="Amplitude_MHCI_rank"),
            AnnotationFactory.build_annotation(value=self.amplitude_mhci_affinity_9mer,
                name="Amplitude_MHCI_affinity_9mer"),
            AnnotationFactory.build_annotation(value=self.amplitude_mhci_MB, name="Amplitude_MHCI_multiple_binding")
        ]

    def get_annotations_mhc2(self) -> List[Annotation]:
        return [
            AnnotationFactory.build_annotation(
                value=self.amplitude_mhcii_affinity, name="Amplitude_MHCII_affinity"),
            AnnotationFactory.build_annotation(value=self.amplitude_mhcii_rank, name="Amplitude_MHCII_rank"),
            AnnotationFactory.build_annotation(value=self.amplitude_mhcii_MB,
                                               name="Amplitude_MHCII_multiple_binding")
        ]
