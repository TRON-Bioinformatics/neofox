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

from neofox.model.neoantigen import PredictedEpitope
from neofox.published_features.Tcell_predictor.tcellpredictor_wrapper import (
    TcellPrediction,
)


class TestTCellPredictor(TestCase):
    def setUp(self) -> None:
        self.tcell_predictor = TcellPrediction()

    # TODO: Franzi maybe you can add some more sensible tests here?

    def test_non_existing_gene(self):
        result = TcellPrediction().calculate_tcell_predictor_score(
            gene="BLAH",
            epitope=PredictedEpitope(wild_type_peptide="BLAHBLAH", mutated_peptide="blaaaah", affinity_mutated=5)
        )
        self.assertEqual(None, result)

    def test_empty_gene(self):
        result = TcellPrediction().calculate_tcell_predictor_score(
            gene=None,
            epitope=PredictedEpitope(wild_type_peptide="BLAHBLAH", mutated_peptide="blaaaah", affinity_mutated=5)
        )
        self.assertEqual(None, result)
        result = TcellPrediction().calculate_tcell_predictor_score(
            gene="",
            epitope=PredictedEpitope(wild_type_peptide="BLAHBLAH", mutated_peptide="blaaaah", affinity_mutated=5)
        )
        self.assertEqual(None, result)
        result = TcellPrediction().calculate_tcell_predictor_score(
            gene="   ",
            epitope=PredictedEpitope(wild_type_peptide="BLAHBLAH", mutated_peptide="blaaaah", affinity_mutated=5)
        )
        self.assertEqual(None, result)

    def test_existing_gene_with_too_short_epitope(self):
        result = TcellPrediction().calculate_tcell_predictor_score(
            gene="BRCA2",
            epitope=PredictedEpitope(wild_type_peptide="CCCCCC", mutated_peptide="CCVCCC", affinity_mutated=5)
        )
        self.assertEqual(None, result)

    def test_existing_gene_with_too_long_epitope(self):
        result = TcellPrediction().calculate_tcell_predictor_score(
            gene="BRCA2",
            epitope=PredictedEpitope(wild_type_peptide="CCCCCCCCCC", mutated_peptide="CCCCCVCCCC", affinity_mutated=5)
        )
        self.assertEqual(None, result)

    def test_existing_gene(self):
        result = TcellPrediction().calculate_tcell_predictor_score(
            gene="BRCA2",
            epitope=PredictedEpitope(wild_type_peptide="CCCCCCCCC", mutated_peptide="CCCCVCCCC", affinity_mutated=5)
        )
        self.assertAlmostEqual(0.30162944008956233, float(result))

    def test_rare_aminoacid(self):
        result = TcellPrediction().calculate_tcell_predictor_score(
            gene="BRCA2",
            epitope=PredictedEpitope(wild_type_peptide="CCCCCCCCC", mutated_peptide="CCCCUCCCC", affinity_mutated=5)
        )
        self.assertIsNone(result)
