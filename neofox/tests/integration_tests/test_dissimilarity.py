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
import os
from unittest import TestCase

import neofox.tests.integration_tests.integration_test_tools as integration_test_tools
from neofox.helpers.blastp_runner import BlastpRunner
from neofox.helpers.runner import Runner
from neofox.published_features.dissimilarity_garnish.dissimilaritycalculator import (
    DissimilarityCalculator,
)


class TestDissimilarity(TestCase):
    def setUp(self):
        self.references, self.configuration = integration_test_tools.load_references()
        self.fastafile = integration_test_tools.create_temp_aminoacid_fasta_file()
        self.runner = Runner()
        self.proteome_blastp_runner = BlastpRunner(
            runner=self.runner, configuration=self.configuration,
            database=os.path.join(self.references.proteome_db, "homo_sapiens"))

    def test_dissimilar_sequences(self):
        result = DissimilarityCalculator(proteome_blastp_runner=self.proteome_blastp_runner).calculate_dissimilarity(
            mhc_mutation="tocino", mhc_affinity=600)
        self.assertEqual(1, result)

    def test_similar_sequences(self):
        result = DissimilarityCalculator(proteome_blastp_runner=self.proteome_blastp_runner).calculate_dissimilarity(
            mhc_mutation="DDDDDD", mhc_affinity=600)
        self.assertTrue(result < 0.000001)

    def test_missing_aminoacid_change(self):
        result = DissimilarityCalculator(proteome_blastp_runner=self.proteome_blastp_runner).calculate_dissimilarity(
            mhc_mutation="DDUDDD", mhc_affinity=600)
        self.assertIsNone(result)

    def test_dissimilarity_mhcii(self):
        # peptide with point mutation
        result = DissimilarityCalculator(proteome_blastp_runner=self.proteome_blastp_runner).calculate_dissimilarity(
            mhc_mutation="LGLSDSQFLQTFLFM", mhc_affinity=430)
        self.assertEqual(result, 0)
        # unsimmilar peptide
        result = DissimilarityCalculator(proteome_blastp_runner=self.proteome_blastp_runner).calculate_dissimilarity(
            mhc_mutation="LELERVLVQY", mhc_affinity=430)
        self.assertEqual(0.0038214427855995936, result)

    def test_dissimilar_sequences(self):
        result = DissimilarityCalculator(proteome_blastp_runner=self.proteome_blastp_runner).calculate_dissimilarity(
            mhc_mutation="tocino", mhc_affinity=600)
        self.assertEqual(1, result)

    def test_affinity_threshold(self):
        # peptide with point mutation
        dissimilariyty_calculator = DissimilarityCalculator(
            proteome_blastp_runner=self.proteome_blastp_runner,
            affinity_threshold=1000
        )
        result = dissimilariyty_calculator.calculate_dissimilarity(
            mhc_mutation="LGLSDSQFLQTFLFM", mhc_affinity=1030, filter_binder=True)
        self.assertIsNone(result)
        result = dissimilariyty_calculator.calculate_dissimilarity(
            mhc_mutation="LGLSDSQFLQTFLFM", mhc_affinity=1030, filter_binder=False)
        self.assertIsNotNone(result)
        result = dissimilariyty_calculator.calculate_dissimilarity(
            mhc_mutation="LGLSDSQFLQTFLFM", mhc_affinity=530, filter_binder=True)
        self.assertIsNotNone(result)

