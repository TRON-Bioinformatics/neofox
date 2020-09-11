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

from neofox.potential_features.aa_index.aa_index import AminoacidIndex


class TestAaIndex(TestCase):

    @classmethod
    def setUpClass(cls):
        cls.aaindex = AminoacidIndex()

    def test_get_aaindex1_frequency(self):
        self.assertEqual(0.946, self.aaindex.aaindex1["KARP850102"]['A'])

    def test_get_aaindex2_frequency(self):
        self.assertEqual(5.7, self.aaindex.aaindex2["VOGG950101"]["A"]["C"])
