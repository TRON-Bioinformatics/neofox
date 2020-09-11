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

from neofox.potential_features.nmer_frequency.nmer_frequency import AminoacidFrequency, FourmerFrequency


class TestNmerFrequency(TestCase):

    @classmethod
    def setUpClass(cls):
        cls.aa_frequencies = AminoacidFrequency()
        cls.fourmer_frequencies = FourmerFrequency()

    def test_get_aminoacid_frequency(self):
        self.assertEqual(7.27339103101, self.aa_frequencies._get_frequency("A"))

    def test_non_existing_aminoacid_frequency(self):
        self.assertIsNone(self.aa_frequencies._get_frequency("X"))

    def test_get_4mer_frequency(self):
        self.assertEqual(0.0148452054899, self.fourmer_frequencies._get_frequency_4mer("AAAAAAAA"))
        self.assertEqual(0.000233496235874, self.fourmer_frequencies._get_frequency_4mer("AAAACCAA"))

    def test_non_existing_4mer_frequency(self):
        self.assertIsNone(self.fourmer_frequencies._get_frequency_4mer("AA"))
        self.assertIsNone(self.fourmer_frequencies._get_frequency_4mer(""))
        self.assertIsNone(self.fourmer_frequencies._get_frequency_4mer("XXXXXXXX"))

    def test_get_4mer_product_frequency(self):
        self.assertEqual(2798.6445730350238, self.aa_frequencies._get_product_4mer_frequencies("AAAAAAAA"))
        self.assertEqual(215.14851896633044, self.aa_frequencies._get_product_4mer_frequencies("AAAACCAA"))

    def test_non_existing_4mer_product_frequency(self):
        self.assertIsNone(self.aa_frequencies._get_product_4mer_frequencies("AA"))
        self.assertIsNone(self.aa_frequencies._get_product_4mer_frequencies(""))
        self.assertIsNone(self.aa_frequencies._get_product_4mer_frequencies("XXXXXXXX"))
