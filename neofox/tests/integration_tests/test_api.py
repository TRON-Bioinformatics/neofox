from typing import List
from unittest import TestCase
from neofox.model.factories import NeoantigenFactory, PatientFactory, NeoepitopeFactory
from neofox.model.neoantigen import Neoantigen
from neofox.neofox import NeoFox
from neofox.references.references import ORGANISM_HOMO_SAPIENS
from neofox.tests.integration_tests import integration_test_tools


class TestApi(TestCase):

    def setUp(self):
        self.references, self.configuration = integration_test_tools.load_references()
        self.hla_database = self.references.get_mhc_database()

    def test_no_mhc2_no_wt(self):

        neoantigen = NeoantigenFactory.build_neoantigen(
            mutated_xmer="SLYNTVATLYCVHQRIDV",
            patient_identifier="123")

        patient = PatientFactory.build_patient(
            identifier="123",
            mhc_alleles=["HLA-A*01:01"],
            mhc_database=self.hla_database)

        annotated_neoantigens = NeoFox(neoantigens=[neoantigen], patients=[patient], num_cpus=2).get_annotations()
        self.assertIsInstance(annotated_neoantigens, List)
        self.assertEqual(len(annotated_neoantigens), 1)
        self.assertIsInstance(annotated_neoantigens[0], Neoantigen)
        self.assertGreater(len(annotated_neoantigens[0].neofox_annotations.annotations), 0)

    def test_no_mhc2(self):

        neoantigen = NeoantigenFactory.build_neoantigen(
            mutated_xmer="SLYNTVATLYCVHQRIDV",
            wild_type_xmer="SLYNTVATAYCVHQRIDV",
            patient_identifier="123")

        patient = PatientFactory.build_patient(
            identifier="123",
            mhc_alleles=["HLA-A*01:01"],
            mhc_database=self.hla_database)

        annotated_neoantigens = NeoFox(neoantigens=[neoantigen], patients=[patient], num_cpus=2).get_annotations()
        self.assertIsInstance(annotated_neoantigens, List)
        self.assertEqual(len(annotated_neoantigens), 1)
        self.assertIsInstance(annotated_neoantigens[0], Neoantigen)
        self.assertGreater(len(annotated_neoantigens[0].neofox_annotations.annotations), 0)

    def test_build_neoepitope_mhc_i(self):

        neoepitope = NeoepitopeFactory.build_neoepitope(
            mutated_peptide="AAAAFAAAA",
            wild_type_peptide="AAAALAAAA",
            allele_mhc_i='HLA-A*01:01',
            organism=ORGANISM_HOMO_SAPIENS,
            mhc_database=self.hla_database
        )
        self.assertIsNotNone(neoepitope)

    def test_build_neoepitope_mhc_i_i(self):

        neoepitope = NeoepitopeFactory.build_neoepitope(
            mutated_peptide="AAAAFAAAA",
            wild_type_peptide="AAAALAAAA",
            isoform_mhc_i_i='DRB1*01:01',
            organism=ORGANISM_HOMO_SAPIENS,
            mhc_database=self.hla_database
        )
        self.assertIsNotNone(neoepitope)

    def test_build_neoepitope_without_mhc(self):

        neoepitope = NeoepitopeFactory.build_neoepitope(
            mutated_peptide="AAAAFAAAA",
            wild_type_peptide="AAAALAAAA",
            patient_identifier='123',
            mhc_database=self.hla_database
        )
        self.assertIsNotNone(neoepitope)
