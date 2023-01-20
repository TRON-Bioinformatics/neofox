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
from logzero import logger
from neofox.published_features.hex.hex import Hex
from neofox.helpers.runner import Runner

import neofox.tests.integration_tests.integration_test_tools as integration_test_tools
from neofox.published_features.hex.pyhex import PyHex


class TestHex(TestCase):
    def setUp(self):
        self.references, self.configuration = integration_test_tools.load_references()
        self.runner = Runner()

    def test_hex(self):
        res = Hex(references=self.references).apply_hex(mut_peptide="FGLAIDVDD")
        self.assertEqual(int(res), 148)

    def test_pyhex(self):
        pyhex = PyHex(iedb_fasta=self.references.get_iedb_fasta())
        res = pyhex.run("FGLAIDVDD")
        self.assertEqual(res, 148)

    def test_comparison(self):
        for i in range(10):
            for k in range(9, 30):
                peptide = integration_test_tools.get_random_kmer(k=k)
                logger.info(peptide)
                res = Hex(references=self.references).apply_hex(mut_peptide=peptide)
                pyhex = PyHex(iedb_fasta=self.references.get_iedb_fasta())
                res_pyhex = pyhex.run(peptide)
                self.assertEqual(float(res), res_pyhex, "Peptide: {}".format(peptide))



