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

import neofox.tests.integration_tests.integration_test_tools as integration_test_tools
from neofox.MHC_predictors.MixMHCpred.abstract_mixmhcpred import AbstractMixMHCpred
from neofox.MHC_predictors.MixMHCpred.mixmhc2pred import MixMhc2Pred
from neofox.MHC_predictors.MixMHCpred.mixmhcpred import MixMHCpred
from neofox.helpers.runner import Runner
from neofox.tests import TEST_HLAI_ALLELES, TEST_HLAII_ALLELES


class TestMixMHCPred(TestCase):

    def setUp(self):
        self.references, self.configuration = integration_test_tools.load_references()
        self.runner = Runner()

    def test_mixmhcpred_epitope_iedb(self):
        mixmhcpred = MixMHCpred(runner=self.runner, configuration=self.configuration)
        # this is an epitope from IEDB of length 9
        mutated = 'NLVPMVATV'
        wild_type = 'NLVPIVATV'
        mixmhcpred.run(sequence_wt=wild_type, sequence_mut=mutated, alleles=TEST_HLAI_ALLELES[0:5])
        self.assertIsNotNone(mixmhcpred.best_peptide)
        self.assertIsNotNone(mixmhcpred.best_score)
        self.assertIsNotNone(mixmhcpred.best_rank)
        self.assertIsNotNone(mixmhcpred.best_allele)

    def test_mixmhcpred_too_small_epitope(self):
        mixmhcpred = MixMHCpred(runner=self.runner, configuration=self.configuration)
        mutated = 'NLVP'
        wild_type = 'NLVP'
        mixmhcpred.run(sequence_wt=wild_type, sequence_mut=mutated, alleles=TEST_HLAI_ALLELES)
        self.assertEqual(None, mixmhcpred.best_peptide)

    def test_mixmhcpred2_epitope_iedb(self):
        mixmhcpred = MixMhc2Pred(runner=self.runner, configuration=self.configuration)
        # this is an epitope from IEDB of length 15
        mutated = 'ENPVVHFFKNIVTPR'
        wild_type = 'ENPVVHIFKNIVTPR'
        mixmhcpred.run(sequence_wt=wild_type, sequence_mut=mutated, alleles=TEST_HLAII_ALLELES)
        self.assertIsNotNone(mixmhcpred.best_peptide)
        self.assertIsNotNone(mixmhcpred.best_rank)
        self.assertIsNotNone(mixmhcpred.best_allele)


    def test_mixmhcpred2_too_small_epitope(self):
        mixmhcpred = MixMhc2Pred(runner=self.runner, configuration=self.configuration)
        # this is an epitope from IEDB of length 15
        mutated = 'ENPVVHFF'
        wild_type = 'ENPVVHFF'
        mixmhcpred.run(sequence_wt=wild_type, sequence_mut=mutated, alleles=TEST_HLAII_ALLELES)
        self.assertEqual(None, mixmhcpred.best_peptide)

    def test_generate_nmers(self):
        result = AbstractMixMHCpred.generate_nmers(
            xmer_wt="DDDDDDDDD", xmer_mut="DDDDDVDDD", lengths=[8, 9, 10, 11])
        self.assertIsNotNone(result)
        self.assertEqual(3, len(result))
        self.assertEqual(1, len(list(filter(lambda x: len(x) == 9, result))))
        self.assertEqual(2, len(list(filter(lambda x: len(x) == 8, result))))
        self.assertEqual(0, len(list(filter(lambda x: len(x) == 7, result))))
        self.assertEqual(0, len(list(filter(lambda x: len(x) == 6, result))))
        self.assertEqual(0, len(list(filter(lambda x: len(x) == 5, result))))
        self.assertEqual(0, len(list(filter(lambda x: len(x) == 4, result))))
        # ['DDDDDVDD', 'DDDDVDDD', 'DDDVDDD', 'DDVDDD', 'DVDDD', 'VDDD', 'DDDDDVDDD', 'DDDDVDDD', 'DDDVDDD', 'DDVDDD',
        #  'DVDDD', 'VDDD']
        logger.debug(result)
