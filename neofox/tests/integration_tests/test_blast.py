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
from neofox.helpers.blastp_runner import (
    BlastpRunner,
)
from logzero import logger
from neofox.helpers.runner import Runner
import neofox.tests.integration_tests.integration_test_tools as integration_test_tools
import os


class TestBlast(TestCase):
    def setUp(self):
        self.references, self.configuration, self.fastafile = self._load_references()
        self.blast_calculator = BlastpRunner(
            runner=Runner(), configuration=self.configuration
        )
        self.proteom_db = self.references.proteome_db


    def _load_references(self):
        references, configuration = integration_test_tools.load_references()
        fastafile = integration_test_tools.create_temp_aminoacid_fasta_file()
        return references, configuration, fastafile

    def test_blast(self):
        sequence = "FIAGLIAIV"
        logger.debug(self.proteom_db)
        database = os.path.join(self.proteom_db, "homo_sapiens")
        logger.debug(database)
        res = self.blast_calculator.get_most_similar_wt_epitope(sequence, database)
        logger.info(res)
