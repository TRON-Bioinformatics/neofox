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

from neofox.model.factories import MhcFactory
from neofox.model.mhc_parser import MhcParser
from neofox.model.neoantigen import Mhc2Name
from neofox.helpers.epitope_helper import EpitopeHelper
import neofox.tests.integration_tests.integration_test_tools as integration_test_tools
from neofox.MHC_predictors.MixMHCpred.mixmhc2pred import MixMHC2pred
from neofox.MHC_predictors.MixMHCpred.mixmhcpred import MixMHCpred
from neofox.helpers.runner import Runner
from neofox.annotation_resources.uniprot.uniprot import Uniprot
from neofox.tests.tools import get_neoantigen


class TestMixMHCPred(TestCase):
    def setUp(self):
        self.references, self.configuration = integration_test_tools.load_references()
        self.runner = Runner()
        mhc_parser = MhcParser.get_mhc_parser(self.references.get_mhc_database())
        self.mixmhcpred = MixMHCpred(
            runner=self.runner, configuration=self.configuration, mhc_parser=mhc_parser
        )
        self.mixmhc2pred = MixMHC2pred(
            runner=self.runner, configuration=self.configuration, mhc_parser=mhc_parser
        )
        self.hla_database = self.references.get_mhc_database()
        self.test_mhc_one = integration_test_tools.get_hla_one_test(self.hla_database)
        self.test_mhc_two = integration_test_tools.get_hla_two_test(self.hla_database)
        self.uniprot = Uniprot(self.references.uniprot_pickle)

    def test_mixmhcpred_epitope_iedb(self):
        # this is an epitope from IEDB of length 9
        mutation = get_neoantigen(mutated_xmer="NLVPMVATV", wild_type_xmer="NLVPIVATV")
        self.mixmhcpred.run(neoantigen=mutation, mhc=self.test_mhc_one, uniprot=self.uniprot)
        best_result = EpitopeHelper.select_best_by_affinity(
            predictions=self.mixmhcpred.results, maximum=True)
        self.assertEquals("NLVPMVATV", best_result.mutated_peptide)
        self.assertAlmostEqual(0.306957, best_result.affinity_mutated, delta=0.00001)
        self.assertEquals(0.6, best_result.rank_mutated)
        self.assertEquals("HLA-A*02:01", best_result.allele_mhc_i.name)

    def test_mixmhcpred_too_small_epitope(self):
        mutation = get_neoantigen(mutated_xmer="NLVP", wild_type_xmer="NLNP")
        self.mixmhcpred.run(neoantigen=mutation, mhc=self.test_mhc_one, uniprot=self.uniprot)
        best_result = EpitopeHelper.select_best_by_affinity(
            predictions=self.mixmhcpred.results, maximum=True)
        self.assertIsNone(best_result.mutated_peptide)
        self.assertIsNone(best_result.rank_mutated)
        self.assertIsNone(best_result.allele_mhc_i.name)
        self.assertIsNone(best_result.affinity_mutated)

    def test_mixmhcpred_not_supported_allele(self):
        """
        this is a combination of neoepitope and HLA alleles from Balachandran
        """
        mutation = get_neoantigen(mutated_xmer="SIYGGLVLI", wild_type_xmer="PIYGGLVLI")
        self.mixmhcpred.run(
            neoantigen=mutation,
            mhc=MhcFactory.build_mhc1_alleles(["A02:01", "B44:02", "C05:17", "C05:01"], self.hla_database),
            uniprot=self.uniprot
        )
        best_result = EpitopeHelper.select_best_by_affinity(
            predictions=self.mixmhcpred.results, maximum=True)
        self.assertEqual('SIYGGLVLI', best_result.mutated_peptide)
        self.assertAlmostEqual(0.158294, best_result.affinity_mutated, places=5)
        self.assertEqual(1, best_result.rank_mutated)
        self.assertEqual('HLA-A*02:01', best_result.allele_mhc_i.name)

    def test_mixmhcpred_rare_aminoacid(self):
        for wild_type_xmer, mutated_xmer in integration_test_tools.mutations_with_rare_aminoacids:
            mutation = get_neoantigen(mutated_xmer=mutated_xmer, wild_type_xmer=wild_type_xmer)
            self.mixmhcpred.run(neoantigen=mutation, mhc=self.test_mhc_one, uniprot=self.uniprot)
            best_result = EpitopeHelper.select_best_by_affinity(
                predictions=self.mixmhcpred.results, maximum=True)
            # rare aminoacids only return empty results when in the mutated sequence
            if EpitopeHelper.contains_rare_amino_acid(mutated_xmer):
                self.assertIsNone(best_result.mutated_peptide)
                self.assertIsNone(best_result.rank_mutated)
                self.assertIsNone(best_result.allele_mhc_i.name)
                self.assertIsNone(best_result.affinity_mutated)
            else:
                self.assertIsNotNone(best_result.mutated_peptide)
                self.assertIsNotNone(best_result.rank_mutated)
                self.assertIsNotNone(best_result.allele_mhc_i)
                self.assertIsNotNone(best_result.affinity_mutated)

    def test_mixmhcpred2_epitope_iedb(self):
        # this is an epitope from IEDB of length 15
        neoantigen = get_neoantigen(
            mutated_xmer="DEVLGEPSQDILVTDQTRLEATISPET",
            wild_type_xmer="DEVLGEPSQDILVIDQTRLEATISPET")
        self.mixmhc2pred.run(
            neoantigen=neoantigen, mhc=self.test_mhc_two,
            uniprot=self.uniprot
        )
        best_result = EpitopeHelper.select_best_by_rank(predictions=self.mixmhc2pred.results)
        self.assertEquals("DEVLGEPSQDILVT", best_result.mutated_peptide)
        self.assertEquals(3.06, best_result.rank_mutated)
        self.assertEquals("HLA-DPA1*01:03-DPB1*04:01", best_result.isoform_mhc_i_i.name)

    def test_mixmhcpred2_epitope_iedb_forcing_no_drb1(self):
        # this is an epitope from IEDB of length 15
        neoantigen = get_neoantigen(
            mutated_xmer="DEVLGEPSQDILVTDQTRLEATISPET",
            wild_type_xmer="DEVLGEPSQDILVIDQTRLEATISPET")
        self.mixmhc2pred.run(
            # forces no DRB1 allele to get as a result one of the composite isoforms
            neoantigen=neoantigen, mhc=[m for m in self.test_mhc_two if m.name != Mhc2Name.DR],
            uniprot=self.uniprot
        )
        best_result = EpitopeHelper.select_best_by_rank(predictions=self.mixmhc2pred.results)
        self.assertEquals("DEVLGEPSQDILVT", best_result.mutated_peptide)
        self.assertEquals(3.06, best_result.rank_mutated)
        self.assertEquals("HLA-DPA1*01:03-DPB1*04:01", best_result.isoform_mhc_i_i.name)

    def test_mixmhcpred2_too_small_epitope(self):
        neoantigen = get_neoantigen(mutated_xmer="ENPVVHFF", wild_type_xmer="ENPVVHFF")
        self.mixmhc2pred.run(neoantigen=neoantigen, mhc=self.test_mhc_two, uniprot=self.uniprot)
        best_result = EpitopeHelper.select_best_by_rank(predictions=self.mixmhc2pred.results)
        self.assertIsNone(best_result.mutated_peptide)
        self.assertIsNone(best_result.rank_mutated)
        self.assertIsNone(best_result.isoform_mhc_i_i.name)

    def test_mixmhcpred2_no_mutation(self):
        for wild_type_xmer, mutated_xmer in integration_test_tools.mutations_with_rare_aminoacids:
            neoantigen = get_neoantigen(mutated_xmer=mutated_xmer, wild_type_xmer=wild_type_xmer)
            self.mixmhc2pred.run(neoantigen=neoantigen, mhc=self.test_mhc_two, uniprot=self.uniprot)
            best_result = EpitopeHelper.select_best_by_rank(predictions=self.mixmhc2pred.results)
            self.assertIsNone(best_result.mutated_peptide)
            self.assertIsNone(best_result.rank_mutated)
            self.assertIsNone(best_result.isoform_mhc_i_i.name)

    def test_mixmhc2pred_rare_aminoacid(self):
        # this is an epitope from IEDB of length 9
        neoantigen = get_neoantigen(mutated_xmer="XTTDSWGKF", wild_type_xmer="XTTDSDGKF")
        self.mixmhc2pred.run(neoantigen=neoantigen, mhc=self.test_mhc_one, uniprot=self.uniprot)
        best_result = EpitopeHelper.select_best_by_rank(predictions=self.mixmhc2pred.results)
        self.assertIsNone(best_result.mutated_peptide)
        self.assertIsNone(best_result.rank_mutated)
        self.assertIsNone(best_result.isoform_mhc_i_i.name)

    def test_mixmhc2pred_allele(self):

        neoantigen = get_neoantigen(mutated_xmer="TNENLDLQELVEKLEKN", wild_type_xmer="TNENLDLQNLVEKLEKN")
        # this is a MHC II genotype which results in no available alleles for MixMHC2pred
        MHC_TWO_NEW = MhcFactory.build_mhc2_alleles(
            [
                "HLA-DRB1*14:54",
                "HLA-DRB1*14:54",
                "HLA-DQA1*01:04",
                "HLA-DQA1*01:04",
                "HLA-DQB1*05:03",
                "HLA-DQB1*05:03",
                "HLA-DPB1*02:01",
                "HLA-DPB1*02:01"
            ],
            self.hla_database
        )
        alleles = self.mixmhc2pred.transform_hla_ii_alleles_for_prediction(MHC_TWO_NEW)
        logger.info(alleles)
        self.mixmhc2pred.run(neoantigen=neoantigen, mhc=MHC_TWO_NEW, uniprot=self.uniprot)
        best_result = EpitopeHelper.select_best_by_rank(predictions=self.mixmhc2pred.results)
        self.assertIsNone(best_result.mutated_peptide)
        self.assertIsNone(best_result.rank_mutated)
        self.assertIsNone(best_result.isoform_mhc_i_i.name)

    def test_generate_nmers(self):
        neoantigen = get_neoantigen(mutated_xmer="DDDDDVDDD", wild_type_xmer="DDDDDDDDD")
        result = EpitopeHelper.generate_nmers(neoantigen=neoantigen, lengths=[8, 9, 10, 11], uniprot=self.uniprot)
        logger.info(result)
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
