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
# Python version with Biopython
# . /etc/profile.d/modules.sh; module load software/python/python-2.7.9

"""
This script takes as neofox table from iCAM pipeline and calculates Literature feature of neoantigens
"""

# import modules
import math
from typing import List

from neofox.model.neoantigen import Annotation, PredictedEpitope
from neofox.model.factories import AnnotationFactory


class VaxRank:

    def logistic_epitope_score(
        self, ic50, midpoint=350.0, width=150.0, ic50_cutoff=5000.0
    ):  # TODO: add these default values into CLI as arguments
        """
        Map from IC50 values to score where 1.0 = strong binder, 0.0 = weak binder
        Default midpoint and width for logistic determined by max likelihood fit
        for data from Alessandro Sette's 1994 paper:
           "The relationship between class I binding affinity
            and immunogenicity of potential cytotoxic T cell epitopes.
        adapted from: https://github.com/openvax/vaxrank/blob/master/vaxrank/epitope_prediction.py
        """
        if ic50 >= ic50_cutoff:
            return 0.0

        rescaled = (float(ic50) - midpoint) / width
        # simplification of 1.0 - logistic(x) = logistic(-x)
        logistic = 1.0 / (1.0 + math.exp(rescaled))

        # since we're scoring IC50 values, let's normalize the output
        # so IC50 near 0.0 always returns a score of 1.0
        normalizer = 1.0 / (1.0 + math.exp(-midpoint / width))

        return logistic / normalizer

    def total_binding(self, epitope_predictions: List[PredictedEpitope]):
        """
        adapted from: https://github.com/openvax/vaxrank/blob/master/vaxrank/epitope_prediction.py
        sums up MHC binding scores of all possible neoepitope candidates, transformed with logistic function into values between 0 and 1
        """
        mut_scores_logistic = 0

        # logistic transformation and sum over all epitopes deriving from mutations
        for p in epitope_predictions:
            mut_scores_logistic += self.logistic_epitope_score(ic50=float(p.affinity_mutated))

        return mut_scores_logistic

    def combined_score(self, expression_score, total_binding_score):
        """
        adapted from: https://github.com/openvax/vaxrank/blob/master/vaxrank/epitope_prediction.py
        final ranking score implemented in VaxRank
        """
        combined_score = None
        try:
            combined_score = float(expression_score) * total_binding_score
        except (ValueError, TypeError):
            pass
        return combined_score

    def get_annotations(self, epitope_predictions: List[PredictedEpitope], expression_score) -> List[Annotation]:
        expression_score = expression_score
        total_binding_score = self.total_binding(epitope_predictions)
        ranking_score = self.combined_score(expression_score=expression_score, total_binding_score=total_binding_score)
        return [
            AnnotationFactory.build_annotation(
                value=total_binding_score, name="Vaxrank_bindingScore"
            ),
            AnnotationFactory.build_annotation(
                value=ranking_score, name="Vaxrank_totalScore"
            ),
        ]
