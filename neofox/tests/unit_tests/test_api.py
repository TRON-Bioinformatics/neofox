from unittest import TestCase

from neofox.model.factories import NeoantigenFactory, PatientFactory
from neofox.model.mhc_parser import MhcParser
from neofox.model.neoantigen import Neoantigen, Patient
from neofox.neofox import NeoFox
from neofox.tests.fake_classes import FakeHlaDatabase
from neofox.tests.integration_tests import integration_test_tools


class TestApi(TestCase):

    def setUp(self) -> None:
        self.hla_database = FakeHlaDatabase()

    def test_no_mhc2_no_wt(self):
        neoantigen = NeoantigenFactory.build_neoantigen(
            mutated_xmer="SLYNTVATLYCVHQRIDV",
            patient_identifier="123")
        self.assertIsInstance(neoantigen, Neoantigen)

        patient = PatientFactory.build_patient(
            identifier="123",
            mhc_alleles=["HLA-A*01:01"],
            mhc_database=self.hla_database)
        self.assertIsInstance(patient, Patient)

    def test_no_mhc2(self):
        neoantigen = NeoantigenFactory.build_neoantigen(
            mutated_xmer="SLYNTVATLYCVHQRIDV",
            wild_type_xmer="SLYNTVATAYCVHQRIDV",
            patient_identifier="123")
        self.assertIsInstance(neoantigen, Neoantigen)

        patient = PatientFactory.build_patient(
            identifier="123",
            mhc_alleles=["HLA-A*01:01"],
            mhc_database=self.hla_database)
        self.assertIsInstance(patient, Patient)

    def test_no_mhc_no_wt(self):
        neoantigen = NeoantigenFactory.build_neoantigen(
            mutated_xmer="SLYNTVATLYCVHQRIDV",
            patient_identifier="123")
        self.assertIsInstance(neoantigen, Neoantigen)

        patient = PatientFactory.build_patient(
            identifier="123",
            mhc2_alleles=["HLA-DRB1*01:01"],
            mhc_database=self.hla_database)
        self.assertIsInstance(patient, Patient)

    def test_no_mhc(self):
        neoantigen = NeoantigenFactory.build_neoantigen(
            mutated_xmer="SLYNTVATLYCVHQRIDV",
            wild_type_xmer="SLYNTVATLACVHQRIDV",
            patient_identifier="123")
        self.assertIsInstance(neoantigen, Neoantigen)

        patient = PatientFactory.build_patient(
            identifier="123",
            mhc2_alleles=["HLA-DRB1*01:01"],
            mhc_database=self.hla_database)
        self.assertIsInstance(patient, Patient)
