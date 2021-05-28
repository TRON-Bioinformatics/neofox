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

from neofox.helpers.blastp_runner import BlastpRunner
from neofox.published_features.neoantigen_fitness.neoantigen_fitness import (
    NeoantigenFitnessCalculator,
)
from neofox.helpers.runner import Runner
import neofox.tests.integration_tests.integration_test_tools as integration_test_tools
from neofox.references.references import IEDB_BLAST_PREFIX


class TestNeoantigenFitness(TestCase):
    def setUp(self):
        self.references, self.configuration, self.fastafile = self._load_references()
        self.iedb_blastp_runner = BlastpRunner(
            runner=Runner(verbose=False), configuration=self.configuration,
            database=os.path.join(self.references.iedb, IEDB_BLAST_PREFIX))
        self.neoantigen_fitness_calculator = NeoantigenFitnessCalculator(iedb_blastp_runner=self.iedb_blastp_runner)

    def _load_references(self):
        references, configuration = integration_test_tools.load_references()
        fastafile = integration_test_tools.create_temp_aminoacid_fasta_file()
        return references, configuration, fastafile

    def test_pathogen_similarity(self):
        # tests a pathogen sequence and expects 1.0 similarity
        result = self.neoantigen_fitness_calculator.get_pathogen_similarity(
            mutation="FIAGLIAIV"
        )
        self.assertEqual(1.0, result)
        # tests a modified pathogen sequence and expects something between 0 and 1
        result = self.neoantigen_fitness_calculator.get_pathogen_similarity(
            mutation="FIAGDAAIV"
        )
        self.assertEqual(0.9923850668756105, result)
        # tests a non pathogen sequence and expects 0 similarity
        result = self.neoantigen_fitness_calculator.get_pathogen_similarity(
            mutation="DDDDDMMDD"
        )
        self.assertEqual(0, result)

    def test_pathogen_similarity_non_existing_aminoacid(self):
        # tests something with a rare aminoacid
        result = self.neoantigen_fitness_calculator.get_pathogen_similarity(
            mutation="DDDDDUDDDD"
        )
        self.assertEqual(0, result)

    def test_amplitude_mhc(self):
        self.assertEqual(
            1.0,
            self.neoantigen_fitness_calculator.calculate_amplitude_mhc(
                score_mutation=1.0, score_wild_type=1.0
            ),
        )
        self.assertEqual(
            0.9997000899730081,
            self.neoantigen_fitness_calculator.calculate_amplitude_mhc(
                score_mutation=1.0, score_wild_type=1.0, apply_correction=True
            ),
        )

    def test_recognition_potential(self):
        self.assertEqual(
            1.0,
            self.neoantigen_fitness_calculator.calculate_recognition_potential(
                amplitude=1.0, pathogen_similarity=1.0, mutation_in_anchor=False
            ),
        )
        self.assertEqual(
            None,
            self.neoantigen_fitness_calculator.calculate_recognition_potential(
                amplitude=None, pathogen_similarity=1.0, mutation_in_anchor=False
            ),
        )

    def test_affinity_threshold(self):
        # tests a pathogen sequence and expects 1.0 similarity
        neoantigen_fitness_calculator = NeoantigenFitnessCalculator(
            iedb_blastp_runner=self.iedb_blastp_runner, affinity_threshold=1
        )
        self.assertIsNone(
            neoantigen_fitness_calculator.calculate_recognition_potential(
                amplitude=1.0, pathogen_similarity=1.0, mhc_affinity_mut=10, mutation_in_anchor=False
            ),
        )
