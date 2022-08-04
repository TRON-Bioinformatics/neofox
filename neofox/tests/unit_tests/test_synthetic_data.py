from unittest import TestCase, skip

import pkg_resources
from faker import Faker

import neofox
from neofox.model.neoantigen import Zygosity, Neoantigen
from neofox.tests.fake_classes import FakeHlaDatabase
from neofox.tests.synthetic_data.factories import PatientProvider, NeoantigenProvider


class TestFactories(TestCase):

    def setUp(self) -> None:
        faker = Faker()
        self.patient_provider = PatientProvider(faker,
            ["HLA-A*74:18", "HLA-A*01:141", "HLA-A*01:12", "HLA-B*07:02", "HLA-B*07:05", "HLA-B*07:06",
             "HLA-C*01:02", "HLA-C*02:10", "HLA-C*03:03"],
            ["DRB1*01:01", "DRB1*01:02", "DRB1*01:03", "DRB1*01:04", "HLA-DPA1*01:01-DPB1*01:01",
             "HLA-DPA1*01:02-DPB1*01:02", "HLA-DPA1*01:03-DPB1*01:03",
             "HLA-DPA1*01:04-DPB1*01:04"
             "HLA-DQA1*01:01-DQB1*01:01", "HLA-DQA1*01:02-DQB1*01:02",
             "HLA-DQA1*01:03-DQB1*01:03", "HLA-DQA1*01:04-DQB1*01:04"],
            FakeHlaDatabase()
                                                )
        self.neoantigen_provider = NeoantigenProvider(faker, proteome_fasta=pkg_resources.resource_filename(
            neofox.tests.__name__, "resources/proteome_test.fa"))

    def test_patient(self):
        patient1 = self.patient_provider.patient()
        patient2 = self.patient_provider.patient()
        self.assertNotEqual(patient1.identifier, patient2.identifier)
        self._assert_patient(patient1)
        self._assert_patient(patient2)

    def _assert_patient(self, patient):
        self.assertIsNotNone(patient.tumor_type)
        for mhc1 in patient.mhc1:
            self.assertTrue(mhc1.zygosity == Zygosity.HETEROZYGOUS)
        for mhc2 in patient.mhc2:
            for g in mhc2.genes:
                self.assertTrue(g.zygosity == Zygosity.HETEROZYGOUS)
        print(patient.to_dict())

    def test_neoantigen(self):
        neoantigen = self.neoantigen_provider.neoantigen()
        self._assert_neoantigen(neoantigen)

    def _assert_neoantigen(self, neoantigen: Neoantigen):
        self.assertIsNotNone(neoantigen.patient_identifier)
        self.assertTrue(len(neoantigen.mutated_xmer) == 27)
        self.assertTrue(len(neoantigen.wild_type_xmer) == 27)
        self.assertNotEqual(neoantigen.wild_type_xmer, neoantigen.mutated_xmer)
        self.assertEqual(neoantigen.wild_type_xmer[0:13], neoantigen.mutated_xmer[0:13])
        self.assertEqual(neoantigen.wild_type_xmer[14:], neoantigen.mutated_xmer[14:])

    def test_neoantigen_and_patient(self):
        patient = self.patient_provider.patient()
        neoantigen = self.neoantigen_provider.neoantigen(patient_identifier=patient.identifier)
        self._assert_patient(patient)
        self._assert_neoantigen(neoantigen)
        self.assertEqual(patient.identifier, neoantigen.patient_identifier)
