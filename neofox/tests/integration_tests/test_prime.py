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
from neofox.helpers.epitope_helper import EpitopeHelper
from neofox.model.factories import MhcFactory
from neofox.model.mhc_parser import HlaParser
import neofox.tests.integration_tests.integration_test_tools as integration_test_tools
from neofox.MHC_predictors.prime import Prime
from neofox.helpers.runner import Runner
from neofox.annotation_resources.uniprot.uniprot import Uniprot
from neofox.tests.tools import get_neoantigen


class TestPrime(TestCase):
    def setUp(self):
        self.references, self.configuration = integration_test_tools.load_references()
        self.runner = Runner()
        self.prime = Prime(
            runner=self.runner, configuration=self.configuration,
            mhc_parser=HlaParser(self.references.get_mhc_database())
        )
        self.hla_database = self.references.get_mhc_database()
        self.test_mhc_one = integration_test_tools.get_hla_one_test(self.hla_database)
        self.uniprot = Uniprot(self.references.uniprot_pickle)

    def test_prime_epitope(self):
        neoantigen = get_neoantigen(mutated_xmer="LVTDQTRLE", wild_type_xmer="LVTDQTRNE")
        self.prime.run(neoantigen=neoantigen, mhc=self.test_mhc_one, uniprot=self.uniprot)
        best_result = EpitopeHelper.select_best_by_affinity(
            predictions=self.prime.results, maximum=True)
        self.assertEquals("LVTDQTRL", best_result.mutated_peptide)
        self.assertAlmostEqual(0.163810, best_result.affinity_mutated, delta=0.00001)
        self.assertEquals(3.00, best_result.rank_mutated)
        self.assertEquals("HLA-C*05:01", best_result.allele_mhc_i.name)

    def test_prime_too_small_epitope(self):
        neoantigen = get_neoantigen(mutated_xmer="NLVP", wild_type_xmer="NLNP")
        self.prime.run(neoantigen=neoantigen, mhc=self.test_mhc_one, uniprot=self.uniprot)
        best_result = EpitopeHelper.select_best_by_affinity(
            predictions=self.prime.results, maximum=True)
        self.assertIsNone(best_result.mutated_peptide)
        self.assertIsNone(best_result.affinity_mutated)
        self.assertIsNone(best_result.rank_mutated)
        self.assertIsNone(best_result.allele_mhc_i.name)

    def test_prime_not_supported_allele(self):
        """
        this is a combination of neoepitope and HLA alleles from Balachandran
        """
        neoantigen = get_neoantigen(mutated_xmer="SIYGGLVLI", wild_type_xmer="PIYGGLVLI")
        self.prime.run(
            neoantigen=neoantigen,
            mhc=MhcFactory.build_mhc1_alleles(["A02:01", "B44:02", "C05:17", "C05:01"], self.hla_database),
            uniprot=self.uniprot
        )
        best_result = EpitopeHelper.select_best_by_affinity(
            predictions=self.prime.results, maximum=True)
        self.assertEqual('SIYGGLVLI', best_result.mutated_peptide)
        self.assertEqual(0.186328, best_result.affinity_mutated)
        self.assertEqual(0.2, best_result.rank_mutated)
        self.assertEqual('HLA-A*02:01', best_result.allele_mhc_i.name)

    def test_prime_rare_aminoacid(self):
        for wild_type_xmer, mutated_xmer in integration_test_tools.mutations_with_rare_aminoacids:
            neoantigen = get_neoantigen(mutated_xmer=mutated_xmer, wild_type_xmer=wild_type_xmer)
            self.prime.run(neoantigen=neoantigen, mhc=self.test_mhc_one, uniprot=self.uniprot)
            # rare aminoacids only return empty results when in the mutated sequence
            best_result = EpitopeHelper.select_best_by_affinity(
                predictions=self.prime.results, maximum=True)
            if EpitopeHelper.contains_rare_amino_acid(mutated_xmer):
                self.assertIsNone(best_result.mutated_peptide)
                self.assertIsNone(best_result.rank_mutated)
                self.assertIsNone(best_result.allele_mhc_i.name)
                self.assertIsNone(best_result.affinity_mutated)
            else:
                self.assertIsNotNone(best_result.mutated_peptide)
                self.assertIsNotNone(best_result.rank_mutated)
                self.assertIsNotNone(best_result.allele_mhc_i.name)
                self.assertIsNotNone(best_result.affinity_mutated)
