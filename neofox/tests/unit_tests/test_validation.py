from unittest import TestCase

from neofox.exceptions import NeofoxDataValidationException
from neofox.model.neoantigen import Gene, Neoantigen, Patient
from neofox.model.validation import ModelValidator
from neofox.tests.fake_classes import FakeAvailableAlleles


class TestModelValidator(TestCase):

    def test_bad_type_raises_exception(self):
        gene = Gene(
            gene="BRCA2",
            transcript_identifier=12345,        # this should be a string instead of an integer
            assembly="hg19")
        self.assertRaises(NeofoxDataValidationException, ModelValidator.validate, gene)

        neoantigen = Neoantigen(
            patient_identifier="1234",
            rna_expression="0.45")              # this should be a float
        self.assertRaises(NeofoxDataValidationException, ModelValidator.validate, neoantigen)

        patient = Patient(
            identifier="1234",
            is_rna_available="Richtig")            # this should be a boolean
        self.assertRaises(NeofoxDataValidationException, ModelValidator.validate, patient)

        # TODO: make validation capture this data types errors!
        gene = Gene(
            gene="BRCA2",
            transcript_identifier=["12345"],    # this should be a string instead of a list of strings
            assembly="hg19")
        ModelValidator.validate(gene)

    def test_good_data_does_not_raise_exceptions(self):
        gene = Gene(
            gene="BRCA2",
            transcript_identifier="12345",
            assembly="hg19")
        ModelValidator.validate(gene)

        neoantigen = Neoantigen(
            patient_identifier="1234",
            rna_expression=0.45)
        ModelValidator.validate(neoantigen)

        patient = Patient(
            identifier="1234",
            is_rna_available=True)
        ModelValidator.validate(patient)

    def test_mhc_i_allele_validation(self):
        patient = Patient(identifier="123", mhc_i_alleles=["HLA-A01:01"])
        valid_patient = ModelValidator.validate_patient(patient)
        self.assertEqual(patient.mhc_i_alleles[0], valid_patient.mhc_i_alleles[0])
        # adds the colon to homogenise representation
        patient.mhc_i_alleles = ["HLA-A0101"]
        valid_patient = ModelValidator.validate_patient(patient)
        self.assertEqual("HLA-A01:01", valid_patient.mhc_i_alleles[0])
        # removes the star
        patient.mhc_i_alleles = ["HLA-A*01:01"]
        valid_patient = ModelValidator.validate_patient(patient)
        self.assertEqual("HLA-A01:01", valid_patient.mhc_i_alleles[0])
        # removes further information
        patient.mhc_i_alleles = ["HLA-A01:01:02:03N"]
        valid_patient = ModelValidator.validate_patient(patient)
        self.assertEqual("HLA-A01:01", valid_patient.mhc_i_alleles[0])
        patient.mhc_i_alleles = ["HLA-A01:01:02N"]
        valid_patient = ModelValidator.validate_patient(patient)
        self.assertEqual("HLA-A01:01", valid_patient.mhc_i_alleles[0])
        patient.mhc_i_alleles = ["HLA-A01:01N"]
        valid_patient = ModelValidator.validate_patient(patient)
        self.assertEqual("HLA-A01:01", valid_patient.mhc_i_alleles[0])

    def test_invalid_mhc_i_alleles(self):
        # P gene is not valid
        patient = Patient(identifier="123", mhc_i_alleles=["HLA-P01:01"])
        self.assertRaises(NeofoxDataValidationException, ModelValidator.validate_patient, patient)
        # serotype 1 is not valid
        patient = Patient(identifier="123", mhc_i_alleles=["HLA-A1:01"])
        self.assertRaises(NeofoxDataValidationException, ModelValidator.validate_patient, patient)
        # no protein information
        patient = Patient(identifier="123", mhc_i_alleles=["HLA-A01"])
        self.assertRaises(NeofoxDataValidationException, ModelValidator.validate_patient, patient)
        # bad protein format
        patient = Patient(identifier="123", mhc_i_alleles=["HLA-A01:ABC"])
        self.assertRaises(NeofoxDataValidationException, ModelValidator.validate_patient, patient)
        # wrong organism, only human
        patient = Patient(identifier="123", mhc_i_alleles=["GOGO-A01:01"])
        self.assertRaises(NeofoxDataValidationException, ModelValidator.validate_patient, patient)
        # no gene
        patient = Patient(identifier="123", mhc_i_alleles=["HLA-0123456"])
        self.assertRaises(NeofoxDataValidationException, ModelValidator.validate_patient, patient)
        # nonsense
        patient = Patient(identifier="123", mhc_i_alleles=["nonsense"])
        self.assertRaises(NeofoxDataValidationException, ModelValidator.validate_patient, patient)

    def test_not_available_mhc_i_alleles(self):
        patient = Patient(identifier="123", mhc_i_alleles=["HLA-A02:01", "HLA-A03:01"])
        valid_patient = ModelValidator.validate_patient(
            patient, available_alleles=FakeAvailableAlleles(available_mch_i=["HLA-A01:01", "HLA-A02:01"]))
        self.assertTrue("HLA-A02:01" in valid_patient.mhc_i_alleles)
        self.assertTrue("HLA-A03:01" not in valid_patient.mhc_i_alleles)
