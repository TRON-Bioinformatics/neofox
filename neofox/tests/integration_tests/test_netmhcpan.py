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
from neofox.helpers import intermediate_files
from neofox.helpers.runner import Runner
from neofox.MHC_predictors.netmhcpan.netmhcIIpan_prediction import NetMhcIIPanPredictor
from neofox.MHC_predictors.netmhcpan.netmhcpan_prediction import NetMhcPanPredictor
from neofox.tests import TEST_MHC_ONE, TEST_MHC_TWO


class TestNetMhcPanPredictor(TestCase):
    def setUp(self):
        references, self.configuration = integration_test_tools.load_references()
        self.runner = Runner()
        self.available_alleles = references.get_available_alleles()

    def test_netmhcpan_epitope_iedb(self):
        netmhcpan_predictor = NetMhcPanPredictor(
            runner=self.runner, configuration=self.configuration
        )
        # this is an epitope from IEDB of length 9
        mutated = "NLVPMVATV"
        predictions = netmhcpan_predictor.mhc_prediction(
            sequence=mutated,
            mhc_alleles=TEST_MHC_ONE,
            set_available_mhc=self.available_alleles.get_available_mhc_i(),
        )
        self.assertEqual(18, len(predictions))

    def test_netmhcpan_too_small_epitope(self):
        netmhcpan_predictor = NetMhcPanPredictor(
            runner=self.runner, configuration=self.configuration
        )
        mutated = "NLVP"
        predictions = netmhcpan_predictor.mhc_prediction(
            sequence=mutated,
            mhc_alleles=TEST_MHC_ONE,
            set_available_mhc=self.available_alleles.get_available_mhc_i(),
        )
        self.assertEqual(0, len(predictions))

    def test_netmhc2pan_epitope_iedb(self):
        netmhc2pan_predictor = NetMhcIIPanPredictor(
            runner=self.runner, configuration=self.configuration
        )
        # this is an epitope from IEDB of length 15
        mutated = "ENPVVHFFKNIVTPR"
        predictions = netmhc2pan_predictor.mhcII_prediction(
            sequence=mutated,
            mhc_alleles=netmhc2pan_predictor.generate_mhc2_alelle_combinations(
                TEST_MHC_TWO
            ),
        )
        self.assertEqual(10, len(predictions))
        for p in predictions:
            self.assertIsNotNone(p.peptide)
            self.assertIsNotNone(p.hla)
            self.assertIsNotNone(p.affinity_score)
            self.assertIsNotNone(p.pos)
            self.assertIsNotNone(p.rank)

    def test_netmhc2pan_too_small_epitope(self):
        netmhc2pan_predictor = NetMhcIIPanPredictor(
            runner=self.runner, configuration=self.configuration
        )
        # this is an epitope from IEDB of length 15
        mutated = "ENPVVH"
        predictions = netmhc2pan_predictor.mhcII_prediction(
            sequence=mutated,
            mhc_alleles=netmhc2pan_predictor.generate_mhc2_alelle_combinations(
                TEST_MHC_TWO
            ),
        )
        self.assertEqual(0, len(predictions))
