from unittest import TestCase
from logzero import logger

import input.tests.integration_tests.integration_test_tools as integration_test_tools
from input.predictors.MixMHCpred.abstract_mixmhcpred import AbstractMixMHCpred
from input.predictors.MixMHCpred.mixmhc2pred import MixMhc2Pred
from input.predictors.MixMHCpred.mixmhcpred import MixMHCpred
from input.helpers.runner import Runner
from input.tests import TEST_HLAI_ALLELES, TEST_HLAII_ALLELES


class TestMixMHCPred(TestCase):

    def setUp(self):
        self.references, self.configuration = integration_test_tools.load_references()
        self.runner = Runner()

    def test_mixmhcpred_epitope_iedb(self):
        mixmhcpred = MixMHCpred(runner=self.runner, configuration=self.configuration)
        # this is an epitope from IEDB of length 9
        mutated = 'NLVPMVATV'
        wild_type = 'NLVPIVATV'
        mixmhcpred.main(xmer_wt=wild_type, xmer_mut=mutated, alleles=TEST_HLAI_ALLELES[0:5])
        self.assertIsNotNone(mixmhcpred.all_peptides)
        self.assertEqual(set("NLVPMVATV|LVPMVATV|NLVPMVAT".split("|")), set(mixmhcpred.all_peptides.split("|")))
        logger.debug(mixmhcpred.all_peptides)
        self.assertIsNotNone(mixmhcpred.all_scores)
        logger.debug(mixmhcpred.all_scores)
        self.assertIsNotNone(mixmhcpred.all_ranks)
        logger.debug(mixmhcpred.all_ranks)
        self.assertEqual(len(mixmhcpred.all_ranks.split('|')), len(mixmhcpred.all_peptides.split('|')))
        self.assertEqual(set("0.306957|-0.479377|-0.522931".split("|")), set(mixmhcpred.all_scores.split("|")))
        self.assertTrue("77" in mixmhcpred.all_ranks.split('|'))
        self.assertIsNotNone(mixmhcpred.all_alleles)
        self.assertIsNotNone(mixmhcpred.best_peptide)
        self.assertIsNotNone(mixmhcpred.best_score)
        self.assertIsNotNone(mixmhcpred.best_rank)
        self.assertIsNotNone(mixmhcpred.best_allele)
        self.assertIsNotNone(mixmhcpred.best_peptide_wt)
        self.assertIsNotNone(mixmhcpred.best_score_wt)
        self.assertIsNotNone(mixmhcpred.best_rank_wt)
        self.assertIsNotNone(mixmhcpred.difference_score_mut_wt)

    def test_mixmhcpred_too_small_epitope(self):
        mixmhcpred = MixMHCpred(runner=self.runner, configuration=self.configuration)
        mutated = 'NLVP'
        wild_type = 'NLVP'
        mixmhcpred.main(xmer_wt=wild_type, xmer_mut=mutated, alleles=TEST_HLAI_ALLELES)
        self.assertEqual("NA", mixmhcpred.all_peptides)

    def test_mixmhcpred2_epitope_iedb(self):
        mixmhcpred = MixMhc2Pred(runner=self.runner, configuration=self.configuration)
        # this is an epitope from IEDB of length 15
        mutated = 'ENPVVHFFKNIVTPR'
        wild_type = 'ENPVVHIFKNIVTPR'
        mixmhcpred.main(xmer_wt=wild_type, xmer_mut=mutated, alleles=TEST_HLAII_ALLELES)
        self.assertIsNotNone(mixmhcpred.all_peptides)
        self.assertTrue("ENPVVHFFKNIVTP" in mixmhcpred.all_peptides.split('|'))
        self.assertTrue("NPVVHFFKNIVTP" in mixmhcpred.all_peptides.split('|'))
        logger.debug(mixmhcpred.all_peptides)
        self.assertIsNotNone(mixmhcpred.all_ranks)
        self.assertEqual(len(mixmhcpred.all_ranks.split('|')), len(mixmhcpred.all_peptides.split('|')))
        self.assertTrue("0.116547" in mixmhcpred.all_ranks.split('|'))
        self.assertTrue("0.276218", mixmhcpred.all_ranks.split('|'))
        logger.debug(mixmhcpred.all_ranks)
        self.assertIsNotNone(mixmhcpred.all_alleles)
        self.assertIsNotNone(mixmhcpred.best_peptide)
        self.assertIsNotNone(mixmhcpred.best_rank)
        self.assertIsNotNone(mixmhcpred.best_allele)
        self.assertIsNotNone(mixmhcpred.best_peptide_wt)
        self.assertIsNotNone(mixmhcpred.best_score_wt)
        self.assertIsNotNone(mixmhcpred.best_rank_wt)
        self.assertIsNotNone(mixmhcpred.difference_score_mut_wt)

    def test_mixmhcpred2_too_small_epitope(self):
        mixmhcpred = MixMhc2Pred(runner=self.runner, configuration=self.configuration)
        # this is an epitope from IEDB of length 15
        mutated = 'ENPVVHFF'
        wild_type = 'ENPVVHFF'
        mixmhcpred.main(xmer_wt=wild_type, xmer_mut=mutated, alleles=TEST_HLAII_ALLELES)
        self.assertEqual("NA", mixmhcpred.all_peptides)

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
