from unittest import TestCase, skip

import pkg_resources
from faker import Faker

import neofox
from neofox.model.neoantigen import Zygosity
from neofox.tests.fake_classes import FakeHlaDatabase
from neofox.tests.synthetic_data.factories import PatientProvider, NeoantigenProvider


class TestFactories(TestCase):

    def setUp(self) -> None:
        faker = Faker()
        self.patient_provider = PatientProvider(faker,
            ["HLA-A*01:01", "HLA-A*01:02", "HLA-A*01:03", "HLA-B*01:01", "HLA-B*01:02", "HLA-B*01:03",
             "HLA-C*01:01", "HLA-C*01:02", "HLA-C*01:03"],
            ["DRB10101", "DRB10102", "DRB10103", "DRB10104", "HLA-DPA10101-DPB10101", "HLA-DPA10102-DPB10102", "HLA-DPA10103-DPB10103",
             "HLA-DPA10104-DPB10104"
             "HLA-DQA10101-DQB10101", "HLA-DQA10102-DQB10102", "HLA-DQA10103-DQB10103", "HLA-DQA10104-DQB10104"],
            FakeHlaDatabase()
                                                )
        self.neoantigen_provider = NeoantigenProvider(faker, proteome_fasta=pkg_resources.resource_filename(
            neofox.tests.__name__, "resources/proteome_test.fa"))

    @skip
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

    def _assert_neoantigen(self, neoantigen):
        self.assertIsNotNone(neoantigen.identifier)
        self.assertIsNotNone(neoantigen.patient_identifier)
        self.assertTrue(len(neoantigen.mutation.mutated_xmer) == 27)
        self.assertTrue(len(neoantigen.mutation.wild_type_xmer) == 27)
        self.assertNotEqual(neoantigen.mutation.wild_type_xmer, neoantigen.mutation.mutated_xmer)
        self.assertEqual(neoantigen.mutation.wild_type_xmer[0:13], neoantigen.mutation.mutated_xmer[0:13])
        self.assertEqual(neoantigen.mutation.wild_type_xmer[14:], neoantigen.mutation.mutated_xmer[14:])

    @skip
    def test_neoantigen_and_patient(self):
        patient = self.patient_provider.patient()
        neoantigen = self.neoantigen_provider.neoantigen(patient_identifier=patient.identifier)
        self._assert_patient(patient)
        self._assert_neoantigen(neoantigen)
        self.assertEqual(patient.identifier, neoantigen.patient_identifier)
