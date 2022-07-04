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
import math
import os
from neofox.model.validation import ModelValidator
from neofox.model.neoantigen import Annotation, PredictedEpitope
from neofox.model.factories import AnnotationFactory


THRESHOLD_IMPROVED_BINDER = 1.2

BETA = 0.11387
BLOSUM62_FILE_NAME = "BLOSUM62-2.matrix.txt"


class SelfSimilarityCalculator:
    def __init__(self):
        blosum_file = os.path.join(
            os.path.abspath(os.path.dirname(__file__)), BLOSUM62_FILE_NAME
        )
        blosum_dict = self._load_blosum(blosum_file)
        self.k1 = self._compute_k1(blosum_dict)

    def _compute_k1(self, blosum_dict):
        K1 = {}
        for i in list(blosum_dict.keys()):
            x = K1.get(i, {})
            for j in list(blosum_dict[i].keys()):
                x[j] = math.pow(blosum_dict[i][j], BETA)
            K1[i] = x
        return K1

    def _load_blosum(self, blosum):
        blosum_dict = {}
        colid = []
        rowid = []
        c = 0
        with open(blosum) as f:
            for line in f:
                c += 1
                if c == 1:
                    colid = line.strip("\n").split(" ")
                    continue
                w = line.strip("\n").split(" ")
                id = w[0]
                v = [float(x) for x in w[1:]]
                rowid.append(id)
                x = blosum_dict.get(id, {})
                for i, vi in enumerate(v):
                    x[colid[i]] = vi
                blosum_dict[id] = x
        return blosum_dict

    def compute_k_hat_3(self, x, y):  # K^3
        return self._compute_k3(x, y) / math.sqrt(
            self._compute_k3(x, x) * self._compute_k3(y, y)
        )

    def _compute_k3(self, f, g):
        max_k = min(len(f), len(g))
        s = 0
        for k in range(1, max_k + 1):
            for i in range(len(f) - (k - 1)):
                u = f[i : i + k]
                for j in range(len(g) - (k - 1)):
                    v = g[j : j + k]
                    s += self._compute_k2k(u, v, self.k1)
        return s

    def _compute_k2k(self, u, v, K1):
        if len(u) != len(v):
            return None
        k = len(u)
        p = K1[u[0]][v[0]]
        for i in range(1, k):
            p = p * K1[u[i]][v[i]]
        return p

    def get_self_similarity(self, mutated_peptide, wt_peptide):
        """
        Returns self-similiarity between mutated and wt epitope according to Bjerregard et al.,
        Argument mhc indicates if determination for MHC I or MHC II epitopes
        """
        self_similarity = None
        if not ModelValidator.has_peptide_rare_amino_acids(mutated_peptide) and \
                not ModelValidator.has_peptide_rare_amino_acids(wt_peptide):
            try:
                self_similarity = str(self.compute_k_hat_3(mutated_peptide, wt_peptide))
            except ZeroDivisionError:
                pass
        return self_similarity

    def is_improved_binder(self, score_mutation, score_wild_type) -> bool:
        """
        This function checks if mutated epitope is improved binder according to Bjerregard et al.
        """
        improved_binder = None
        try:
            improved_binder = (
                score_wild_type / score_mutation >= THRESHOLD_IMPROVED_BINDER
            )
        except (ZeroDivisionError, ValueError, TypeError):
            pass
        return improved_binder

    def self_similarity_of_conserved_binder_only(
        self, is_improved_binder, similarity
    ):
        """
        this function returns selfsimilarity for conserved binder but not for improved binder
        """
        result = None
        try:
            if not is_improved_binder:
                result = similarity
        except (ZeroDivisionError, ValueError):
            pass
        return result

    def get_annnotations(
            self, epitope_mhci: PredictedEpitope, epitope_mhcii: PredictedEpitope) -> List[Annotation]:

        improved_binding_mhci = None
        self_similarity_mhci = None
        self_similarity_mhcii = None
        if epitope_mhci and epitope_mhci.mutated_peptide and epitope_mhci.wild_type_peptide:
            improved_binding_mhci = self.is_improved_binder(
                score_mutation=epitope_mhci.rank_mutated,
                score_wild_type=epitope_mhci.rank_wild_type,
            )
            self_similarity_mhci = self.get_self_similarity(
                mutated_peptide=epitope_mhci.mutated_peptide,
                wt_peptide=epitope_mhci.wild_type_peptide,
            )
        if epitope_mhcii and epitope_mhcii.mutated_peptide and epitope_mhcii.wild_type_peptide:
            self_similarity_mhcii = self.get_self_similarity(
                mutated_peptide=epitope_mhcii.mutated_peptide,
                wt_peptide=epitope_mhcii.wild_type_peptide,
            )
        annotations = [
            AnnotationFactory.build_annotation(
                value=improved_binding_mhci, name="Improved_Binder_MHCI"
            ),
            AnnotationFactory.build_annotation(
                value=self_similarity_mhcii, name="Selfsimilarity_MHCII"
            ),
            AnnotationFactory.build_annotation(
                value=self_similarity_mhci, name="Selfsimilarity_MHCI"
            ),
            AnnotationFactory.build_annotation(
                value=self.self_similarity_of_conserved_binder_only(
                    is_improved_binder=improved_binding_mhci,
                    similarity=self_similarity_mhci,
                ),
                name="Selfsimilarity_MHCI_conserved_binder",
            ),
        ]
        return annotations

    def get_annotations_epitope_mhcii(self, epitope: PredictedEpitope) -> List[Annotation]:
        return [
            AnnotationFactory.build_annotation(
                value=self.get_self_similarity(
                    mutated_peptide=epitope.mutated_peptide, wt_peptide=epitope.wild_type_peptide),
                name='Selfsimilarity')
            ]

    def get_annotations_epitope_mhci(self, epitope: PredictedEpitope) -> List[Annotation]:
        is_improved_binder = self.is_improved_binder(
            score_mutation=epitope.rank_mutated, score_wild_type=epitope.rank_wild_type)
        self_similarity = self.get_self_similarity(
            mutated_peptide=epitope.mutated_peptide, wt_peptide=epitope.wild_type_peptide)

        return [
            AnnotationFactory.build_annotation(value=is_improved_binder, name='Improved_Binder_MHCI'),
            AnnotationFactory.build_annotation(value=self_similarity, name='Selfsimilarity'),
            AnnotationFactory.build_annotation(
                value=self.self_similarity_of_conserved_binder_only(
                    similarity=self_similarity, is_improved_binder=is_improved_binder),
                name='Selfsimilarity_conserved_binder')
        ]