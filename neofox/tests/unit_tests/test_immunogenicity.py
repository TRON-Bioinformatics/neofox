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

from neofox.published_features.iedb_immunogenicity.iedb import IEDBimmunogenicity


class TestImmunogenicity(TestCase):
    def setUp(self):
        self.immunogenicity_calculator = IEDBimmunogenicity()

    def test_immunogenicity(self):
        result = self.immunogenicity_calculator.calculate_iedb_immunogenicity(
            epitope="ENPVVHFF", mhc_allele="HLA-A*68:01", mhc_score=600
        )
        self.assertGreater(result, 0)
        result = self.immunogenicity_calculator.calculate_iedb_immunogenicity(
            epitope="ENPVVHFF",
            mhc_allele="HLA-A*68:01",
            mhc_score=600,
            affin_filtering=True,
        )
        self.assertIsNone(result)
