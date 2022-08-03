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
from neofox.helpers.blastp_runner import BlastpRunner
from neofox.helpers.runner import Runner
from neofox.MHC_predictors.netmhcpan.netmhcIIpan_prediction import NetMhcIIPanPredictor
from neofox.MHC_predictors.netmhcpan.netmhcpan_prediction import NetMhcPanPredictor
from neofox.model.mhc_parser import MhcParser
from neofox.references.references import ORGANISM_MUS_MUSCULUS


class TestNetMhcPanPredictorOnMouse(TestCase):
    def setUp(self):
        references, self.configuration = integration_test_tools.load_references(organism=ORGANISM_MUS_MUSCULUS)
        self.runner = Runner()
        self.available_alleles = references.get_available_alleles()
        self.test_mhc_one = integration_test_tools.get_h2_one_test(references.get_mhc_database())
        self.test_mhc_two = integration_test_tools.get_h2_two_test(references.get_mhc_database())
        self.mhc_parser = MhcParser.get_mhc_parser(references.get_mhc_database())
        self.proteome_blastp_runner = BlastpRunner(
            runner=self.runner, configuration=self.configuration,
            database=references.get_proteome_database())

    def test_netmhcpan_epitope_iedb(self):
        netmhcpan_predictor = NetMhcPanPredictor(
            runner=self.runner, configuration=self.configuration, mhc_parser=self.mhc_parser,
            blastp_runner=self.proteome_blastp_runner
        )
        # this is an epitope from IEDB of length 9
        mutated = "NLVPMVATV"
        available_alleles = netmhcpan_predictor.get_only_available_alleles(
            mhc_alleles=self.test_mhc_one, set_available_mhc=self.available_alleles.get_available_mhc_i())
        predictions = netmhcpan_predictor.mhc_prediction(sequence=mutated, available_alleles=available_alleles)
        self.assertEqual(9, len(predictions))

    def test_netmhcpan_too_small_epitope(self):
        netmhcpan_predictor = NetMhcPanPredictor(
            runner=self.runner, configuration=self.configuration, mhc_parser=self.mhc_parser,
            blastp_runner=self.proteome_blastp_runner
        )
        mutated = "NLVP"
        available_alleles = netmhcpan_predictor.get_only_available_alleles(
            mhc_alleles=self.test_mhc_one, set_available_mhc=self.available_alleles.get_available_mhc_i())
        predictions = netmhcpan_predictor.mhc_prediction(sequence=mutated, available_alleles=available_alleles)
        self.assertEqual(0, len(predictions))

    def test_netmhcpan_rare_aminoacid(self):
        netmhcpan_predictor = NetMhcPanPredictor(
            runner=self.runner, configuration=self.configuration, mhc_parser=self.mhc_parser,
            blastp_runner=self.proteome_blastp_runner
        )
        # this is an epitope from IEDB of length 9
        mutated = "XTTDSWGKF"
        available_alleles = netmhcpan_predictor.get_only_available_alleles(
            mhc_alleles=self.test_mhc_one, set_available_mhc=self.available_alleles.get_available_mhc_i())
        predictions = netmhcpan_predictor.mhc_prediction(sequence=mutated, available_alleles=available_alleles)
        self.assertEqual(9, len(predictions))

    def test_netmhc2pan_epitope_iedb(self):
        netmhc2pan_predictor = NetMhcIIPanPredictor(
            runner=self.runner, configuration=self.configuration, mhc_parser=self.mhc_parser,
            blastp_runner=self.proteome_blastp_runner
        )
        # this is an epitope from IEDB of length 15
        mutated = "ENPVVHFFKNIVTPR"
        combinations = netmhc2pan_predictor.represent_mhc2_isoforms(
            netmhc2pan_predictor.generate_mhc2_alelle_combinations(self.test_mhc_two))
        predictions = netmhc2pan_predictor.mhc2_prediction(sequence=mutated, mhc_alleles=combinations)
        self.assertEqual(2, len(predictions))
        self._assert_predictions(predictions)

    def _assert_predictions(self, predictions):
        for p in predictions:
            self.assertIsNotNone(p.mutated_peptide)
            self.assertIsNotNone(p.allele_mhc_i)
            self.assertIsNotNone(p.affinity_mutated)
            self.assertIsNotNone(p.position)
            self.assertIsNotNone(p.rank_mutated)

    def test_netmhc2pan_too_small_epitope(self):
        netmhc2pan_predictor = NetMhcIIPanPredictor(
            runner=self.runner, configuration=self.configuration, mhc_parser=self.mhc_parser,
            blastp_runner=self.proteome_blastp_runner
        )
        # this is an epitope from IEDB of length 15
        mutated = "ENPVVH"
        combinations = netmhc2pan_predictor.represent_mhc2_isoforms(
            netmhc2pan_predictor.generate_mhc2_alelle_combinations(self.test_mhc_two))
        predictions = netmhc2pan_predictor.mhc2_prediction(sequence=mutated, mhc_alleles=combinations)
        self.assertEqual(0, len(predictions))
        self._assert_predictions(predictions)

    def test_netmhc2pan_rare_aminoacid(self):
        netmhc2pan_predictor = NetMhcIIPanPredictor(
            runner=self.runner, configuration=self.configuration, mhc_parser=self.mhc_parser,
            blastp_runner=self.proteome_blastp_runner
        )
        # this is an epitope from IEDB of length 15
        mutated = "XTTDSWGKFDDDDDDDDD"
        combinations = netmhc2pan_predictor.represent_mhc2_isoforms(
            netmhc2pan_predictor.generate_mhc2_alelle_combinations(self.test_mhc_two))
        predictions = netmhc2pan_predictor.mhc2_prediction(sequence=mutated, mhc_alleles=combinations)
        self.assertEqual(8, len(predictions))
        self._assert_predictions(predictions)
