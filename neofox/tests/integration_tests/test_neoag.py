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

from neofox.model.neoantigen import Annotation, PredictedEpitope, MhcAllele
from neofox.published_features.neoag.neoag_gbm_model import NeoagCalculator
from neofox.helpers.runner import Runner
import neofox.tests.integration_tests.integration_test_tools as integration_test_tools
from neofox.tests.tools import get_neoantigen


class TestNeoantigenFitness(TestCase):
    def setUp(self):
        self.references, self.configuration = integration_test_tools.load_references()
        self.fastafile = integration_test_tools.create_temp_aminoacid_fasta_file()
        self.runner = Runner()

    def test_neoag(self):

        mutation = get_neoantigen(
                mutated_xmer=  "DEVLGEPSQDILVTDQTRLEATISPET",
                wild_type_xmer="DEVLGEPSQDILVIDQTRLEATISPET"
        )
        result = NeoagCalculator(
            runner=self.runner, configuration=self.configuration
        ).get_annotation(
            epitope_mhci=PredictedEpitope(
                mutated_peptide="ILVTDQTRL", wild_type_peptide="ILVIDQTRL",
                affinity_mutated=0, position=0, allele_mhc_i=MhcAllele(name="hla"), rank_mutated=0
            ),
            neoantigen=mutation,
        )
        self.assertTrue(isinstance(result, Annotation))
        self.assertTrue(float(result.value) > 0)

    def test_affinity_threshold(self):
        mutation = get_neoantigen(
                mutated_xmer="DEVLGEPSQDILVTDQTRLEATISPET",
                wild_type_xmer="DEVLGEPSQDILVIDQTRLEATISPET",
        )
        result = NeoagCalculator(
            runner=self.runner, configuration=self.configuration
        ).get_annotation(
            epitope_mhci=PredictedEpitope(
                mutated_peptide="DDDDDV", affinity_mutated=10, position=0, allele_mhc_i=MhcAllele(name="hla"),
                rank_mutated=0
            ),
            neoantigen=mutation
        )
        self.assertEqual(result.value, "NA")
