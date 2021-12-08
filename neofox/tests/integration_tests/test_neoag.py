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

from neofox.model.neoantigen import Annotation
from neofox.published_features.neoag.neoag_gbm_model import NeoagCalculator
from neofox.helpers.runner import Runner
import neofox.tests.integration_tests.integration_test_tools as integration_test_tools
from neofox.MHC_predictors.netmhcpan.abstract_netmhcpan_predictor import PredictedEpitope
from neofox.tests.tools import get_mutation


class TestNeoantigenFitness(TestCase):
    def setUp(self):
        self.references, self.configuration = integration_test_tools.load_references()
        self.fastafile = integration_test_tools.create_temp_aminoacid_fasta_file()
        self.runner = Runner()

    def test_neoag(self):

        mutation = get_mutation(
                mutated_xmer="DEVLGEPSQDILVTDQTRLEATISPET",
                wild_type_xmer="DEVLGEPSQDILVIDQTRLEATISPET"
        )
        result = NeoagCalculator(
            runner=self.runner, configuration=self.configuration
        ).get_annotation(
            sample_id="12345",
            mutated_peptide_mhci=PredictedEpitope(
                peptide="DDDDDV", affinity_score=0, pos=0, hla="hla", rank=0
            ),
            wt_peptide_mhci=PredictedEpitope(
                peptide="DDDDDD", affinity_score=0, pos=0, hla="hla", rank=0
            ),
            mutation=mutation,
            peptide_variant_position="123"
        )
        self.assertTrue(isinstance(result, Annotation))
        self.assertTrue(float(result.value) > 0)

    def test_affinity_threshold(self):
        mutation = get_mutation(
                mutated_xmer="DEVLGEPSQDILVTDQTRLEATISPET",
                wild_type_xmer="DEVLGEPSQDILVIDQTRLEATISPET",
        )
        result = NeoagCalculator(
            runner=self.runner, configuration=self.configuration, affinity_threshold=1
        ).get_annotation(
            sample_id="12345",
            mutated_peptide_mhci=PredictedEpitope(
                peptide="DDDDDV", affinity_score=10, pos=0, hla="hla", rank=0
            ),
            wt_peptide_mhci=PredictedEpitope(
                peptide="DDDDDD", affinity_score=0, pos=0, hla="hla", rank=0
            ),
            mutation=mutation,
            peptide_variant_position="123"
        )
        self.assertEqual(result.value, "NA")
