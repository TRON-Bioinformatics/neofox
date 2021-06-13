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

from neofox.helpers.epitope_helper import EpitopeHelper
from neofox.model.conversion import ModelValidator, ModelConverter
from neofox.model.mhc_parser import MhcParser
from neofox.model.neoantigen import Mutation

import neofox.tests.integration_tests.integration_test_tools as integration_test_tools
from neofox.published_features.prime import Prime
from neofox.helpers.runner import Runner
from neofox.annotation_resources.uniprot.uniprot import Uniprot


class TestPrime(TestCase):
    def setUp(self):
        self.references, self.configuration = integration_test_tools.load_references()
        self.runner = Runner()
        self.prime = Prime(
            runner=self.runner, configuration=self.configuration,
            mhc_parser=MhcParser(self.references.get_hla_database())
        )
        self.hla_database = self.references.get_hla_database()
        self.test_mhc_one = integration_test_tools.get_mhc_one_test(self.hla_database)
        self.uniprot = Uniprot(self.references.uniprot_pickle)

    def test_prime_epitope(self):
        mutation = ModelValidator._validate_mutation(
            Mutation(mutated_xmer="LVTDQTRLE", wild_type_xmer="LVTDQTRNE")
        )
        best_peptide, best_rank, best_allele, best_score = self.prime.run(
            mutation=mutation, mhc=self.test_mhc_one, uniprot=self.uniprot
        )
        self.assertEquals("LVTDQTRL", best_peptide)
        self.assertAlmostEqual(0.163810, best_score, delta=0.00001)
        self.assertEquals(3.00, best_rank)
        self.assertEquals("HLA-C*05:01", best_allele)

    def test_prime_too_small_epitope(self):
        mutation = ModelValidator._validate_mutation(
            Mutation(mutated_xmer="NLVP", wild_type_xmer="NLNP")
        )
        best_peptide, best_rank, best_allele, best_score = self.prime.run(
            mutation=mutation, mhc=self.test_mhc_one, uniprot=self.uniprot
        )
        self.assertIsNone(best_peptide)
        self.assertIsNone(best_score)
        self.assertIsNone(best_rank)
        self.assertIsNone(best_allele)

    def test_prime_not_supported_allele(self):
        """
        this is a combination of neoepitope and HLA alleles from Balachandran
        """
        mutation = ModelValidator._validate_mutation(
            Mutation(mutated_xmer="SIYGGLVLI", wild_type_xmer="PIYGGLVLI")
        )
        best_peptide, best_rank, best_allele, best_score = self.prime.run(
            mutation=mutation,
            mhc=ModelConverter.parse_mhc1_alleles(["A02:01", "B44:02", "C05:17", "C05:01"], self.hla_database),
            uniprot=self.uniprot
        )
        self.assertEqual('SIYGGLVLI', best_peptide)
        self.assertEqual(0.186328, best_score)
        self.assertEqual(0.2, best_rank)
        self.assertEqual('HLA-A*02:01', best_allele)

    def test_prime_rare_aminoacid(self):
        for wild_type_xmer, mutated_xmer in integration_test_tools.mutations_with_rare_aminoacids:
            mutation = ModelValidator._validate_mutation(
                Mutation(mutated_xmer=mutated_xmer, wild_type_xmer=wild_type_xmer)
            )
            best_peptide, best_rank, best_allele, best_score = self.prime.run(
                mutation=mutation, mhc=self.test_mhc_one, uniprot=self.uniprot
            )
            # rare aminoacids only return empty results when in the mutated sequence
            if EpitopeHelper.contains_rare_amino_acid(mutated_xmer):
                self.assertIsNone(best_peptide)
                self.assertIsNone(best_rank)
                self.assertIsNone(best_allele)
                self.assertIsNone(best_score)
            else:
                self.assertIsNotNone(best_peptide)
                self.assertIsNotNone(best_rank)
                self.assertIsNotNone(best_allele)
                self.assertIsNotNone(best_score)
