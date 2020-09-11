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
from collections import defaultdict
from unittest import TestCase

from neofox.published_features.neoantigen_fitness.neoantigen_fitness import NeoantigenFitnessCalculator
from neofox.helpers.runner import Runner
import neofox.tests.integration_tests.integration_test_tools as integration_test_tools


class TestNeoantigenFitness(TestCase):

    def setUp(self):
        self.references, self.configuration, self.fastafile = self._load_references()
        self.neoantigen_fitness_calculator = NeoantigenFitnessCalculator(
            runner=Runner(), configuration=self.configuration, iedb=self.references.iedb)

    def _load_references(self):
        references, configuration = integration_test_tools.load_references()
        fastafile = integration_test_tools.create_temp_aminoacid_fasta_file()
        return references, configuration, fastafile

    def test_pathogen_similarity(self):
        # tests a pathogen sequence and expects 1.0 similarity
        result = self.neoantigen_fitness_calculator.wrap_pathogen_similarity(mutation='FIAGLIAIV')
        self.assertEqual(1.0, result)
        # tests a modified pathogen sequence and expects something between 0 and 1
        result = self.neoantigen_fitness_calculator.wrap_pathogen_similarity(mutation='FIAGDAAIV')
        self.assertLess(float(result), 1.0)
        self.assertGreater(float(result), 0.0)
        # tests a non pathogen sequence and expects 0 similarity
        result = self.neoantigen_fitness_calculator.wrap_pathogen_similarity(mutation='DDDDDMMDD')
        self.assertEqual(0, result)

    def test_amplitude_mhc(self):
        self.assertEqual(1.0, self.neoantigen_fitness_calculator.calculate_amplitude_mhc(
            score_mutation=1.0, score_wild_type=1.0))
        self.assertEqual(0.9997000899730081, self.neoantigen_fitness_calculator.calculate_amplitude_mhc(
            score_mutation=1.0, score_wild_type=1.0, apply_correction=True))

    def test_recognition_potential(self):
        self.assertEqual(1.0, self.neoantigen_fitness_calculator.calculate_recognition_potential(
            amplitude=1.0, pathogen_similarity=1.0, mutation_in_anchor=False))
        self.assertEqual(None, self.neoantigen_fitness_calculator.calculate_recognition_potential(
            amplitude=None, pathogen_similarity=1.0, mutation_in_anchor=False))
