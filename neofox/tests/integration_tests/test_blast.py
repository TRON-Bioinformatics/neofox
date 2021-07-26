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
import time
from unittest import TestCase

import numpy as np

from neofox.helpers.blastp_runner import (
    BlastpRunner,
)
from logzero import logger
from neofox.helpers.runner import Runner
import neofox.tests.integration_tests.integration_test_tools as integration_test_tools


class TestBlast(TestCase):
    def setUp(self):
        self.references, self.configuration, self.fastafile = self._load_references()
        self.proteome_blastp_runner = BlastpRunner(
            runner=Runner(verbose=False), configuration=self.configuration,
            database=self.references.get_proteome_database())
        self.iedb_blastp_runner = BlastpRunner(
            runner=Runner(verbose=False), configuration=self.configuration,
            database=self.references.get_iedb_database())

    def _load_references(self):
        references, configuration = integration_test_tools.load_references()
        fastafile = integration_test_tools.create_temp_aminoacid_fasta_file()
        return references, configuration, fastafile

    def test_blast(self):
        sequence = "FIAGLIAIV"
        res = self.proteome_blastp_runner.get_most_similar_wt_epitope(sequence)
        self.assertIsNotNone(res)
        self.assertEqual('FVAGLIVLL', res)

    def test_calculate_similarity_database(self):
        sequence = "FIAGLIAIV"
        res = self.iedb_blastp_runner.calculate_similarity_database(sequence)
        self.assertIsNotNone(res)
        self.assertEqual(1.0, res)

    def test_performance_most_similar_wt_epitope(self):
        times = []
        for k in [9, 15, 25, 50]:
            for _ in range(10):
                kmer = integration_test_tools.get_random_kmer(k=9)
                start = time.time()
                similar_wt = self.proteome_blastp_runner.get_most_similar_wt_epitope(kmer)
                times.append(time.time() - start)
                logger.info("{}: {}".format(kmer, similar_wt))
            logger.info("Average time for {}-mers: {}".format(k, np.mean(times)))
            logger.info("Standard deviation for {}-mers: {}".format(k, np.std(times)))

    def test_performance_calculate_similarity_database(self):
        times = []
        for k in [9, 15, 25, 50]:
            for _ in range(10):
                kmer = integration_test_tools.get_random_kmer(k=9)
                start = time.time()
                similarity = self.proteome_blastp_runner.calculate_similarity_database(kmer)
                times.append(time.time() - start)
                logger.info("{}: {}".format(kmer, similarity))
            logger.info("Average time for {}-mers: {}".format(k, np.mean(times)))
            logger.info("Standard deviation for {}-mers: {}".format(k, np.std(times)))
