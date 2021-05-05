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

from neofox.model.conversion import ModelValidator
from neofox.model.neoantigen import Annotation
from neofox.model.wrappers import AnnotationFactory
from neofox.MHC_predictors.netmhcpan.combine_netmhcpan_pred_multiple_binders import (
    BestAndMultipleBinder,
)

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

    def get_self_similarity(self, mutation, wild_type):
        """
        Returns self-similiarity between mutated and wt epitope according to Bjerregard et al.,
        Argument mhc indicates if determination for MHC I or MHC II epitopes
        """
        self_similarity = None
        if not ModelValidator.has_peptide_rare_amino_acids(mutation) and \
                not ModelValidator.has_peptide_rare_amino_acids(wild_type):
            try:
                self_similarity = str(self.compute_k_hat_3(mutation, wild_type))
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
        self, has_conserved_binder, similarity
    ):
        """
        this function returns selfsimilarity for conserved binder but not for improved binder
        """
        result = None
        try:
            # TODO: is this logic correct? improved binder is synonymous to conserved binder or opposite?
            if not has_conserved_binder:
                result = similarity
        except (ZeroDivisionError, ValueError):
            pass
        return result

    def get_annnotations(self, netmhcpan: BestAndMultipleBinder) -> List[Annotation]:

        improved_binding_mhc1 = None
        self_similarity_mhc1 = None
        if netmhcpan.best_epitope_by_rank.peptide and netmhcpan.best_wt_epitope_by_rank.peptide:
            improved_binding_mhc1 = self.is_improved_binder(
                score_mutation=netmhcpan.best_epitope_by_rank.rank,
                score_wild_type=netmhcpan.best_wt_epitope_by_rank.rank,
            )
            self_similarity_mhc1 = self.get_self_similarity(
                mutation=netmhcpan.best_epitope_by_rank.peptide,
                wild_type=netmhcpan.best_wt_epitope_by_rank.peptide,
            )
        annotations = [
            AnnotationFactory.build_annotation(
                value=improved_binding_mhc1, name="Improved_Binder_MHCI"
            ),
            AnnotationFactory.build_annotation(
                value=self.self_similarity_of_conserved_binder_only(
                    has_conserved_binder=improved_binding_mhc1,
                    similarity=self_similarity_mhc1,
                ),
                name="Selfsimilarity_MHCI_conserved_binder",
            ),
        ]
        return annotations
