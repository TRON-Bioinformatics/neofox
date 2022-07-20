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
from neofox.published_features.priority_score import PriorityScore


class TestPriorityScore(TestCase):
    def setUp(self):
        self.priority_calculator = PriorityScore()

    def test_priority(self):
        result = self.priority_calculator.calc_priority_score(
            vaf_dna=0.35,
            vaf_rna=0.33,
            transcript_expr=12,
            no_mismatch=1,
            score_mut=1.1,
            score_wt=10,
            mut_not_in_prot=True,
        )
        self.assertGreater(result, 0)
        result = self.priority_calculator.calc_priority_score(
            vaf_dna=None,
            vaf_rna=0.33,
            transcript_expr=12,
            no_mismatch=1,
            score_mut=1.1,
            score_wt=10,
            mut_not_in_prot=True,
        )
        self.assertGreater(result, 0)
        result = self.priority_calculator.calc_priority_score(
            vaf_dna=0.35,
            vaf_rna=None,
            transcript_expr=12,
            no_mismatch=1,
            score_mut=1.1,
            score_wt=10,
            mut_not_in_prot=True,
        )
        self.assertGreater(result, 0)
        result = self.priority_calculator.calc_priority_score(
            vaf_dna=None,
            vaf_rna=-1,
            transcript_expr=12,
            no_mismatch=1,
            score_mut=1.1,
            score_wt=10,
            mut_not_in_prot=True,
        )
        self.assertEqual(result, None)
        result = self.priority_calculator.calc_priority_score(
            vaf_dna=0.35,
            vaf_rna=0.33,
            transcript_expr=None,
            no_mismatch=1,
            score_mut=1.1,
            score_wt=10,
            mut_not_in_prot=True,
        )
        self.assertEqual(result, None)
        result = self.priority_calculator.calc_priority_score(
            vaf_dna=0.35,
            vaf_rna=0.33,
            transcript_expr=None,
            no_mismatch=1,
            score_mut=1.1,
            score_wt=10,
            mut_not_in_prot=True,
        )
        self.assertEqual(result, None)
