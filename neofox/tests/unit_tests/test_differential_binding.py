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

from neofox.published_features.differential_binding.differential_binding import (
    DifferentialBinding,
)


class TestDifferentialBinding(TestCase):
    def setUp(self):
        self.diffbdg_calculator = DifferentialBinding(affinity_threshold=500)

    def test_dai(self):
        result = self.diffbdg_calculator.dai(
            score_mutation=50, score_wild_type=500, affin_filtering=False
        )
        self.assertEqual(result, 450)
        result = self.diffbdg_calculator.dai(
            score_mutation=550, score_wild_type=2000, affin_filtering=True
        )
        self.assertEqual(result, None)
        result = self.diffbdg_calculator.dai(
            score_mutation=50, score_wild_type=10, affin_filtering=True
        )
        self.assertLess(result, 0.0)

    def test_classify_adn_cdn(self):
        result = self.diffbdg_calculator.classify_adn_cdn(
            score_mutation=2,
            amplitude=11,
            bdg_cutoff_classical=50,
            bdg_cutoff_alternative=5000,
            amplitude_cutoff=10,
            category="CDN",
        )
        self.assertEqual(result, True)
        result = self.diffbdg_calculator.classify_adn_cdn(
            score_mutation=70,
            amplitude=11,
            bdg_cutoff_classical=50,
            bdg_cutoff_alternative=5000,
            amplitude_cutoff=10,
            category="CDN",
        )
        self.assertEqual(result, False)
        result = self.diffbdg_calculator.classify_adn_cdn(
            score_mutation=70,
            amplitude=11,
            bdg_cutoff_classical=50,
            bdg_cutoff_alternative=5000,
            amplitude_cutoff=10,
            category="ADN",
        )
        self.assertEqual(result, True)
        result = self.diffbdg_calculator.classify_adn_cdn(
            score_mutation=70,
            amplitude=5,
            bdg_cutoff_classical=50,
            bdg_cutoff_alternative=5000,
            amplitude_cutoff=10,
            category="ADN",
        )
        self.assertEqual(result, False)
        result = self.diffbdg_calculator.classify_adn_cdn(
            score_mutation=None,
            amplitude=5,
            bdg_cutoff_classical=50,
            bdg_cutoff_alternative=5000,
            amplitude_cutoff=10,
            category="ADN",
        )
        self.assertEqual(result, None)

    def test_affinity_threshold(self):
        diffbdg_calculator = DifferentialBinding(affinity_threshold=1000)
        result = diffbdg_calculator.dai(
            score_mutation=50, score_wild_type=530, affin_filtering=False
        )
        self.assertIsNotNone(result)
        result = diffbdg_calculator.dai(
            score_mutation=550, score_wild_type=530, affin_filtering=True
        )
        self.assertIsNotNone(result)
        result = diffbdg_calculator.dai(
            score_mutation=1030, score_wild_type=1030, affin_filtering=True
        )
        self.assertIsNone(result)
