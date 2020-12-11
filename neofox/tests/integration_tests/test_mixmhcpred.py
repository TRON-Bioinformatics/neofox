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

from neofox.model.conversion import ModelValidator, ModelConverter
from neofox.model.neoantigen import Mutation

from neofox.helpers.epitope_helper import EpitopeHelper

import neofox.tests.integration_tests.integration_test_tools as integration_test_tools
from neofox.MHC_predictors.MixMHCpred.mixmhc2pred import MixMhc2Pred
from neofox.MHC_predictors.MixMHCpred.mixmhcpred import MixMHCpred
from neofox.helpers.runner import Runner
from neofox.tests import TEST_MHC_ONE, TEST_MHC_TWO


class TestMixMHCPred(TestCase):
    def setUp(self):
        self.references, self.configuration = integration_test_tools.load_references()
        self.runner = Runner()
        self.mixmhcpred = MixMHCpred(
            runner=self.runner, configuration=self.configuration
        )
        self.mixmhc2pred = MixMhc2Pred(
            runner=self.runner, configuration=self.configuration
        )

    def test_mixmhcpred_epitope_iedb(self):
        # this is an epitope from IEDB of length 9
        mutation = ModelValidator._validate_mutation(
            Mutation(mutated_xmer="NLVPMVATV", wild_type_xmer="NLVPIVATV")
        )
        best_peptide, best_rank, best_allele, best_score = self.mixmhcpred.run(
            mutation=mutation, mhc=TEST_MHC_ONE
        )
        self.assertEquals("NLVPMVATV", best_peptide)
        self.assertAlmostEqual(0.306957, best_score, delta=0.00001)
        self.assertEquals(0.6, best_rank)
        self.assertEquals("A0201", best_allele)

    def test_mixmhcpred_too_small_epitope(self):
        mutation = ModelValidator._validate_mutation(
            Mutation(mutated_xmer="NLVP", wild_type_xmer="NLNP")
        )
        best_peptide, best_rank, best_allele, best_score = self.mixmhcpred.run(
            mutation=mutation, mhc=TEST_MHC_ONE
        )
        self.assertIsNone(best_peptide)
        self.assertIsNone(best_score)
        self.assertIsNone(best_rank)
        self.assertIsNone(best_allele)

    def test_mixmhcpred_no_mutation(self):
        mutation = ModelValidator._validate_mutation(
            Mutation(mutated_xmer="NNNNNNNNN", wild_type_xmer="NNNNNNNNN")
        )
        best_peptide, best_rank, best_allele, best_score = self.mixmhcpred.run(
            mutation=mutation, mhc=TEST_MHC_ONE
        )
        self.assertIsNone(best_peptide)
        self.assertIsNone(best_score)
        self.assertIsNone(best_rank)
        self.assertIsNone(best_allele)

    def test_mixmhcpred_not_supported_allele(self):
        """
        this is a combination of neoepitope and HLA alleles from Balachandran
        """
        mutation = ModelValidator._validate_mutation(
            Mutation(mutated_xmer="SIYGGLVLI", wild_type_xmer="PIYGGLVLI")
        )
        best_peptide, best_rank, best_allele, best_score = self.mixmhcpred.run(
            mutation=mutation, mhc=ModelConverter.parse_mhc1_alleles(["A02:01", "B44:02", "C05:17", "C05:01"])
        )
        self.assertEqual('SIYGGLVLI', best_peptide)
        self.assertEqual(0.15829400000000002, best_score)
        self.assertEqual(1, best_rank)
        self.assertEqual('A0201', best_allele)

    def test_mixmhcpred_rare_aminoacid(self):
        # this is an epitope from IEDB of length 9
        mutation = ModelValidator._validate_mutation(
            Mutation(mutated_xmer="XTTDSWGKF", wild_type_xmer="XTTDSDGKF")
        )
        best_peptide, best_rank, best_allele, best_score = self.mixmhcpred.run(
            mutation=mutation, mhc=TEST_MHC_ONE
        )
        self.assertIsNone(best_peptide)
        self.assertIsNone(best_rank)
        self.assertIsNone(best_allele)
        self.assertIsNone(best_score)

    def test_mixmhcpred2_epitope_iedb(self):
        # this is an epitope from IEDB of length 15
        mutation = ModelValidator._validate_mutation(
            Mutation(mutated_xmer="ENPVVHFFKNIVTPR", wild_type_xmer="ENPVVHIFKNIVTPR")
        )
        best_peptide, best_rank, best_allele = self.mixmhc2pred.run(
            mutation=mutation, mhc=TEST_MHC_TWO
        )
        self.assertEquals("NPVVHFFKNIVTPR", best_peptide)
        self.assertEquals(2.16, best_rank)
        self.assertEquals("DRB1_08_01", best_allele)

    def test_mixmhcpred2_too_small_epitope(self):
        mutation = ModelValidator._validate_mutation(
            Mutation(mutated_xmer="ENPVVHFF", wild_type_xmer="ENPVVHFF")
        )
        best_peptide, best_rank, best_allele = self.mixmhc2pred.run(
            mutation=mutation, mhc=TEST_MHC_TWO
        )
        self.assertIsNone(best_peptide)
        self.assertIsNone(best_rank)
        self.assertIsNone(best_allele)

    def test_mixmhcpred2_no_mutation(self):
        # this is an epitope from IEDB of length 15
        mutation = ModelValidator._validate_mutation(
            Mutation(mutated_xmer="ENPVVHFFKNIVTPR", wild_type_xmer="ENPVVHFFKNIVTPR")
        )
        best_peptide, best_rank, best_allele = self.mixmhc2pred.run(
            mutation=mutation, mhc=TEST_MHC_TWO
        )
        self.assertIsNone(best_peptide)
        self.assertIsNone(best_rank)
        self.assertIsNone(best_allele)

    def test_mixmhc2pred_rare_aminoacid(self):
        # this is an epitope from IEDB of length 9
        mutation = ModelValidator._validate_mutation(
            Mutation(mutated_xmer="XTTDSWGKF", wild_type_xmer="XTTDSDGKF")
        )
        best_peptide, best_rank, best_allele = self.mixmhc2pred.run(
            mutation=mutation, mhc=TEST_MHC_ONE
        )
        self.assertIsNone(best_peptide)
        self.assertIsNone(best_rank)
        self.assertIsNone(best_allele)

    def test_generate_nmers(self):
        mutation = ModelValidator._validate_mutation(
            Mutation(mutated_xmer="DDDDDVDDD", wild_type_xmer="DDDDDDDDD")
        )
        result = EpitopeHelper.generate_nmers(mutation=mutation, lengths=[8, 9, 10, 11])
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
