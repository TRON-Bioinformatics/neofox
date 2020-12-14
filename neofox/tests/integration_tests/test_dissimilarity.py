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

import neofox.tests.integration_tests.integration_test_tools as integration_test_tools
from neofox.helpers.runner import Runner
from neofox.published_features.dissimilarity_garnish.dissimilaritycalculator import (
    DissimilarityCalculator,
)


class TestDissimilarity(TestCase):
    def setUp(self):
        self.references, self.configuration = integration_test_tools.load_references()
        self.fastafile = integration_test_tools.create_temp_aminoacid_fasta_file()
        self.runner = Runner()

    def test_dissimilar_sequences(self):
        result = DissimilarityCalculator(
            runner=self.runner,
            configuration=self.configuration,
            proteome_db=self.references.proteome_db,
        ).calculate_dissimilarity(mhc_mutation="tocino", mhc_affinity="velocidad")
        self.assertEqual(1, result)

    def test_similar_sequences(self):
        result = DissimilarityCalculator(
            runner=self.runner,
            configuration=self.configuration,
            proteome_db=self.references.proteome_db,
        ).calculate_dissimilarity(mhc_mutation="DDDDDD", mhc_affinity="DDDDDD")
        self.assertTrue(result < 0.000001)

    def test_missing_aminoacid_change(self):
        result = DissimilarityCalculator(
            runner=self.runner,
            configuration=self.configuration,
            proteome_db=self.references.proteome_db,
        ).calculate_dissimilarity(mhc_mutation="DDUDDD", mhc_affinity="DDYDDD")
        self.assertIsNone(result)

    def test_dissimilarity_mhcii(self):
        # peptide with point mutation
        result = DissimilarityCalculator(
            runner=self.runner,
            configuration=self.configuration,
            proteome_db=self.references.proteome_db,
        ).calculate_dissimilarity(mhc_mutation="LGLSDSQFLQTFLFM", mhc_affinity="430")
        self.assertEqual(result, 0)
        # unsimmilar peptide
        result = DissimilarityCalculator(
            runner=self.runner,
            configuration=self.configuration,
            proteome_db=self.references.proteome_db,
        ).calculate_dissimilarity(mhc_mutation="LFTSPIMTKSAEMIV", mhc_affinity="430")
        self.assertEqual(9.713825893542527e-06, result)
