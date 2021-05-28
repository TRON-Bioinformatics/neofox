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
from logzero import logger
from neofox.helpers.blastp_runner import BlastpRunner
from neofox.model.neoantigen import Annotation
from neofox.model.wrappers import AnnotationFactory
from neofox.MHC_predictors.netmhcpan.combine_netmhcpan_pred_multiple_binders import (
    BestAndMultipleBinder,
)
from neofox import AFFINITY_THRESHOLD_DEFAULT
from neofox.published_features.differential_binding.amplitude import Amplitude


class NeoantigenFitnessCalculator:

    def __init__(self, iedb_blastp_runner: BlastpRunner, affinity_threshold=AFFINITY_THRESHOLD_DEFAULT):
        self.affinity_threshold = affinity_threshold
        self.iedb_blastp_runner = iedb_blastp_runner

    def get_pathogen_similarity(self, mutation):
        pathsim = self.iedb_blastp_runner.calculate_similarity_database(peptide=mutation)
        logger.info(
            "Peptide {} has a pathogen similarity of {}".format(mutation, pathsim)
        )
        return pathsim

    def calculate_amplitude_mhc(
        self, score_mutation, score_wild_type, apply_correction=False
    ):
        """
        This function calculates the amplitude between mutated and wt epitope according to Balachandran et al.
        when affinity is used, use correction from Luksza et al. *1/(1+0.0003*aff_wt)
        """
        amplitude_mhc = None
        try:
            candidate_amplitude_mhc = score_wild_type / score_mutation
            if apply_correction:  # nine_mer or affinity:
                amplitude_mhc = candidate_amplitude_mhc * self._calculate_correction(
                    score_wild_type
                )
            else:
                amplitude_mhc = candidate_amplitude_mhc
        except (ZeroDivisionError, ValueError, TypeError):
            pass
        return amplitude_mhc

    def _calculate_correction(self, score_wild_type):
        return 1 / (1 + 0.0003 * float(score_wild_type))

    def calculate_recognition_potential(
        self,
        amplitude: float,
        pathogen_similarity: float,
        mutation_in_anchor: bool,
        mhc_affinity_mut: float = None,
    ):
        """
        This function calculates the recognition potential, defined by the product of amplitude and pathogensimiliarity
        of an epitope according to Balachandran et al.
        F_alpha = - max (A_i x R_i)

        Returns (A_i x R_i) value only for nonanchor mutation and epitopes of length 9; only considered by Balachandran
        """
        recognition_potential = None
        try:
            candidate_recognition_potential = amplitude * pathogen_similarity
            if mhc_affinity_mut:
                if not mutation_in_anchor and mhc_affinity_mut < self.affinity_threshold:
                    recognition_potential = candidate_recognition_potential
            else:
                if not mutation_in_anchor:
                    recognition_potential = candidate_recognition_potential
        except (ValueError, TypeError):
            pass
        return recognition_potential

    def get_annotations(
        self, netmhcpan: BestAndMultipleBinder, amplitude: Amplitude
    ) -> List[Annotation]:
        pathogen_similarity_9mer = None
        recognition_potential = None
        if netmhcpan.best_ninemer_epitope_by_affinity.peptide:
            pathogen_similarity_9mer = self.get_pathogen_similarity(
                mutation=netmhcpan.best_ninemer_epitope_by_affinity.peptide
            )
            if pathogen_similarity_9mer is not None:
                recognition_potential = self.calculate_recognition_potential(
                            amplitude=amplitude.amplitude_mhci_affinity_9mer,
                            pathogen_similarity=pathogen_similarity_9mer,
                            mutation_in_anchor=netmhcpan.mutation_in_anchor_9mer,
                            mhc_affinity_mut=netmhcpan.best_ninemer_epitope_by_affinity.affinity_score,
                        )
        annotations = [
            AnnotationFactory.build_annotation(
                name="Pathogensimiliarity_MHCI_affinity_9mer",
                value=pathogen_similarity_9mer,
            ),
            AnnotationFactory.build_annotation(
                name="Recognition_Potential_MHCI_affinity_9mer",
                value=recognition_potential
            ),
        ]
        return annotations
