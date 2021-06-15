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
from neofox.published_features.Tcell_predictor.tcellpredictor_wrapper import (
    TcellPrediction,
)


class TestTCellPredictor(TestCase):
    def setUp(self) -> None:
        self.tcell_predictor = TcellPrediction()

    # TODO: Franzis maybe you can add some more sensible tests here?

    def test_non_existing_gene(self):
        result = TcellPrediction(affinity_threshold=10)._calculate_tcell_predictor_score(
            gene="BLAH",
            substitution="blaaaah",
            epitope="BLAHBLAH",
            score=5
        )
        self.assertEqual(None, result)

    def test_empty_gene(self):
        result = TcellPrediction(affinity_threshold=10)._calculate_tcell_predictor_score(
            gene=None,
            substitution="blaaaah",
            epitope="BLAHBLAH",
            score=5
        )
        self.assertEqual(None, result)
        result = TcellPrediction(affinity_threshold=10)._calculate_tcell_predictor_score(
            gene="",
            substitution="blaaaah",
            epitope="BLAHBLAH",
            score=5
        )
        self.assertEqual(None, result)
        result = TcellPrediction(affinity_threshold=10)._calculate_tcell_predictor_score(
            gene="   ",
            substitution="blaaaah",
            epitope="BLAHBLAH",
            score=5
        )
        self.assertEqual(None, result)

    def test_existing_gene_with_too_short_epitope(self):
        result = TcellPrediction(affinity_threshold=10)._calculate_tcell_predictor_score(
            gene="BRCA2", substitution="C", epitope="CCCCCC", score=5
        )
        self.assertEqual(None, result)

    def test_existing_gene_with_too_long_epitope(self):
        result = TcellPrediction(affinity_threshold=10)._calculate_tcell_predictor_score(
            gene="BRCA2", substitution="C", epitope="CCCCCCCCCC", score=5
        )
        self.assertEqual(None, result)

    def test_existing_gene(self):
        result = TcellPrediction(affinity_threshold=10)._calculate_tcell_predictor_score(
            gene="BRCA2",
            substitution="CCCCVCCCC",
            epitope="CCCCCCCCC",
            score=5
        )
        self.assertAlmostEqual(0.2453409331088489, float(result))

    def test_rare_aminoacid(self):
        result = TcellPrediction(affinity_threshold=10)._calculate_tcell_predictor_score(
            gene="BRCA2",
            substitution="CU",
            epitope="CCCCUCCCC",
            score=5
        )
        self.assertIsNone(result)

    def test_affinity_threshold(self):
        result = TcellPrediction(affinity_threshold=1)._calculate_tcell_predictor_score(
            gene="BRCA2",
            substitution="CCCCVCCCC",
            epitope="CCCCCCCCC",
            score=5
        )
        self.assertIsNone(result)
