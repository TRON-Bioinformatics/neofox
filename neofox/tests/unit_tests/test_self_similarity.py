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
from unittest import TestCase

from neofox.published_features.self_similarity.self_similarity import (
    SelfSimilarityCalculator,
)


class TestSelfSimilarity(TestCase):
    def test_get_self_similarity(self):
        result = SelfSimilarityCalculator().get_self_similarity(
            wild_type="DDD", mutation="DDD"
        )
        self.assertEqual("1.0", result)

    def test_is_improved_binder(self):
        result = SelfSimilarityCalculator().is_improved_binder(
            score_mutation=1.0, score_wild_type=1.3
        )
        self.assertTrue(result)

    def test_compute_self_similarity_calculator(self):

        s = SelfSimilarityCalculator()
        self.assertEqual(s.compute_k_hat_3("AAAAA", "AAAAA"), 1.0)
        for i in range(5):
            self.assertTrue(s.compute_k_hat_3("AAAAA", "WWWWW" * (i + 1)) < 1.0)
        for i in list(s.k1.keys()):
            if i == "A":
                self.assertTrue(s.compute_k_hat_3("AAAAA", "AA" + i + "AA") == 1.0)
            else:
                self.assertTrue(s.compute_k_hat_3("AAAAA", "AA" + i + "AA") < 1.0)

    def test_rare_aminoacids(self):
        result = SelfSimilarityCalculator().get_self_similarity(
            wild_type="DDD", mutation="DUD"
        )
        self.assertEqual(None, result)
