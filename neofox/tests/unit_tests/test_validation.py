from unittest import TestCase

from neofox.exceptions import NeofoxDataValidationException
from neofox.model.neoantigen import Gene, Neoantigen, Patient, MhcAllele
from neofox.model.validation import ModelValidator


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
        patient = Patient(identifier="123", mhc_i_alleles=[MhcAllele(name="HLA-A01:01")])
        valid_patient = ModelValidator.validate_patient(patient)
        self.assertEqual(patient.mhc_i_alleles[0], valid_patient.mhc_i_alleles[0])
        # adds the colon to homogenise representation
        patient.mhc_i_alleles = [MhcAllele(name="HLA-A0101")]
        valid_patient = ModelValidator.validate_patient(patient)
        self.assertEqual("HLA-A*01:01", valid_patient.mhc_i_alleles[0].name)
        # removes the star
        patient.mhc_i_alleles = [MhcAllele(name="HLA-A*01:01")]
        valid_patient = ModelValidator.validate_patient(patient)
        self.assertEqual("HLA-A*01:01", valid_patient.mhc_i_alleles[0].name)
        # removes further information
        patient.mhc_i_alleles = [MhcAllele(name="HLA-A01:01:02:03N")]
        valid_patient = ModelValidator.validate_patient(patient)
        self.assertEqual("HLA-A*01:01", valid_patient.mhc_i_alleles[0].name)
        patient.mhc_i_alleles = [MhcAllele(name="HLA-A01:01:02N")]
        valid_patient = ModelValidator.validate_patient(patient)
        self.assertEqual("HLA-A*01:01", valid_patient.mhc_i_alleles[0].name)
        patient.mhc_i_alleles = [MhcAllele(name="HLA-A01:01N")]
        valid_patient = ModelValidator.validate_patient(patient)
        self.assertEqual("HLA-A*01:01", valid_patient.mhc_i_alleles[0].name)
        patient.mhc_i_alleles = [MhcAllele(gene="A", group="01", protein="01")]
        valid_patient = ModelValidator.validate_patient(patient)
        self.assertEqual("HLA-A*01:01", valid_patient.mhc_i_alleles[0].name)

    def test_mhc_ii_allele_validation(self):
        patient = Patient(identifier="123",
                          mhc_i_i_alleles=[MhcAllele(name="HLA-DPA101:01"), MhcAllele(name="HLA-DPB101:01")])
        valid_patient = ModelValidator.validate_patient(patient)
        self.assertEqual(patient.mhc_i_i_alleles[0], valid_patient.mhc_i_i_alleles[0])
        # adds the colon to homogenise representation
        patient.mhc_i_i_alleles = [MhcAllele(name="HLA-DPA10101"), MhcAllele(name="HLA-DPB101:01")]
        valid_patient = ModelValidator.validate_patient(patient)
        self.assertEqual("HLA-DPA1*01:01", valid_patient.mhc_i_i_alleles[0].name)
        # removes the star
        patient.mhc_i_i_alleles = [MhcAllele(name="HLA-DPA1*01:01"), MhcAllele(name="HLA-DPB101:01")]
        valid_patient = ModelValidator.validate_patient(patient)
        self.assertEqual("HLA-DPA1*01:01", valid_patient.mhc_i_i_alleles[0].name)
        # removes further information
        patient.mhc_i_i_alleles = [MhcAllele(name="HLA-DPA101:01:02:03N"), MhcAllele(name="HLA-DPB101:01")]
        valid_patient = ModelValidator.validate_patient(patient)
        self.assertEqual("HLA-DPA1*01:01", valid_patient.mhc_i_i_alleles[0].name)
        patient.mhc_i_i_alleles = [MhcAllele(name="HLA-DPA101:01:02N"), MhcAllele(name="HLA-DPB101:01")]
        valid_patient = ModelValidator.validate_patient(patient)
        self.assertEqual("HLA-DPA1*01:01", valid_patient.mhc_i_i_alleles[0].name)
        patient.mhc_i_i_alleles = [MhcAllele(name="HLA-DPA101:01N"), MhcAllele(name="HLA-DPB101:01")]
        valid_patient = ModelValidator.validate_patient(patient)
        self.assertEqual("HLA-DPA1*01:01", valid_patient.mhc_i_i_alleles[0].name)
        patient.mhc_i_i_alleles = [MhcAllele(gene="DPA1", group="01", protein="01"), MhcAllele(name="HLA-DPB101:01")]
        valid_patient = ModelValidator.validate_patient(patient)
        self.assertEqual("HLA-DPA1*01:01", valid_patient.mhc_i_i_alleles[0].name)

    def test_invalid_mhc_i_alleles(self):
        # P gene is not valid
        patient = Patient(identifier="123", mhc_i_alleles=[MhcAllele(name="HLA-P01:01")])
        self.assertRaises(NeofoxDataValidationException, ModelValidator.validate_patient, patient)
        # serotype 1 is not valid
        patient = Patient(identifier="123", mhc_i_alleles=[MhcAllele(name="HLA-A1:01")])
        self.assertRaises(NeofoxDataValidationException, ModelValidator.validate_patient, patient)
        # no protein information
        patient = Patient(identifier="123", mhc_i_alleles=[MhcAllele(name="HLA-A01")])
        self.assertRaises(NeofoxDataValidationException, ModelValidator.validate_patient, patient)
        # bad protein format
        patient = Patient(identifier="123", mhc_i_alleles=[MhcAllele(name="HLA-A01:ABC")])
        self.assertRaises(NeofoxDataValidationException, ModelValidator.validate_patient, patient)
        # wrong organism, only human
        patient = Patient(identifier="123", mhc_i_alleles=[MhcAllele(name="GOGO-A01:01")])
        self.assertRaises(NeofoxDataValidationException, ModelValidator.validate_patient, patient)
        # no gene
        patient = Patient(identifier="123", mhc_i_alleles=[MhcAllele(name="HLA-0123456")])
        self.assertRaises(NeofoxDataValidationException, ModelValidator.validate_patient, patient)
        # nonsense
        patient = Patient(identifier="123", mhc_i_alleles=[MhcAllele(name="nonsense")])
        self.assertRaises(NeofoxDataValidationException, ModelValidator.validate_patient, patient)
        # missing protein
        patient = Patient(identifier="123", mhc_i_alleles=[MhcAllele(gene="A", group="01")])
        self.assertRaises(NeofoxDataValidationException, ModelValidator.validate_patient, patient)
        # bad protein
        patient = Patient(identifier="123", mhc_i_alleles=[MhcAllele(gene="A", group="01", protein="NaN")])
        self.assertRaises(NeofoxDataValidationException, ModelValidator.validate_patient, patient)
        # bad group
        patient = Patient(identifier="123", mhc_i_alleles=[MhcAllele(gene="A", group="NaN", protein="01")])
        self.assertRaises(NeofoxDataValidationException, ModelValidator.validate_patient, patient)
        # non existing gene
        patient = Patient(identifier="123", mhc_i_alleles=[MhcAllele(gene="Z", group="01", protein="01")])
        self.assertRaises(NeofoxDataValidationException, ModelValidator.validate_patient, patient)
        # MHC I non classical
        patient = Patient(identifier="123", mhc_i_alleles=[MhcAllele(gene="E", group="01", protein="01")])
        self.assertRaises(NeofoxDataValidationException, ModelValidator.validate_patient, patient)
        patient = Patient(identifier="123", mhc_i_alleles=[MhcAllele(gene="F", group="01", protein="01")])
        self.assertRaises(NeofoxDataValidationException, ModelValidator.validate_patient, patient)
        patient = Patient(identifier="123", mhc_i_alleles=[MhcAllele(gene="G", group="01", protein="01")])
        self.assertRaises(NeofoxDataValidationException, ModelValidator.validate_patient, patient)

    def test_invalid_mhc_ii_alleles(self):
        # P gene is not valid
        patient = Patient(identifier="123", mhc_i_i_alleles=[MhcAllele(name="HLA-DPR01:01")])
        self.assertRaises(NeofoxDataValidationException, ModelValidator.validate_patient, patient)
        # serotype 1 is not valid
        patient = Patient(identifier="123", mhc_i_i_alleles=[MhcAllele(name="HLA-DPA11:01")])
        self.assertRaises(NeofoxDataValidationException, ModelValidator.validate_patient, patient)
        # no protein information
        patient = Patient(identifier="123", mhc_i_i_alleles=[MhcAllele(name="HLA-DPA101")])
        self.assertRaises(NeofoxDataValidationException, ModelValidator.validate_patient, patient)
        # bad protein format
        patient = Patient(identifier="123", mhc_i_i_alleles=[MhcAllele(name="HLA-DPA101:ABC")])
        self.assertRaises(NeofoxDataValidationException, ModelValidator.validate_patient, patient)
        # wrong organism, only human
        patient = Patient(identifier="123", mhc_i_i_alleles=[MhcAllele(name="GOGO-DPA101:01")])
        self.assertRaises(NeofoxDataValidationException, ModelValidator.validate_patient, patient)
        # no gene
        patient = Patient(identifier="123", mhc_i_i_alleles=[MhcAllele(name="HLA-0123456")])
        self.assertRaises(NeofoxDataValidationException, ModelValidator.validate_patient, patient)
        # nonsense
        patient = Patient(identifier="123", mhc_i_i_alleles=[MhcAllele(name="nonsense")])
        self.assertRaises(NeofoxDataValidationException, ModelValidator.validate_patient, patient)
        # missing protein
        patient = Patient(identifier="123", mhc_i_i_alleles=[MhcAllele(gene="DPA1", group="01")])
        self.assertRaises(NeofoxDataValidationException, ModelValidator.validate_patient, patient)
        # bad protein
        patient = Patient(identifier="123", mhc_i_i_alleles=[MhcAllele(gene="DPA1", group="01", protein="NaN")])
        self.assertRaises(NeofoxDataValidationException, ModelValidator.validate_patient, patient)
        # bad group
        patient = Patient(identifier="123", mhc_i_i_alleles=[MhcAllele(gene="DPA1", group="NaN", protein="01")])
        self.assertRaises(NeofoxDataValidationException, ModelValidator.validate_patient, patient)
        # non existing gene
        patient = Patient(identifier="123", mhc_i_i_alleles=[MhcAllele(gene="DPA1ZZZZ", group="01", protein="01")])
        self.assertRaises(NeofoxDataValidationException, ModelValidator.validate_patient, patient)

    def test_invalid_mhc_i_genotype(self):
        # 3 alleles for gene A
        patient = Patient(
            identifier="123",
            mhc_i_alleles=[MhcAllele(name="HLA-A01:01"), MhcAllele(name="HLA-A01:02"), MhcAllele(name="HLA-A01:03")])
        self.assertRaises(NeofoxDataValidationException, ModelValidator.validate_patient, patient)

    def test_invalid_mhc_ii_genotype(self):
        # 3 alleles for gene A
        patient = Patient(
            identifier="123",
            mhc_i_i_alleles=[MhcAllele(name="HLA-DPA101:01"), MhcAllele(name="HLA-DPA101:02"),
                             MhcAllele(name="HLA-DPA101:03")])
        self.assertRaises(NeofoxDataValidationException, ModelValidator.validate_patient, patient)

    def test_empty_patient_identifier(self):
        patient = Patient(identifier=None)
        self.assertRaises(NeofoxDataValidationException, ModelValidator.validate_patient, patient)
        patient = Patient(identifier="")
        self.assertRaises(NeofoxDataValidationException, ModelValidator.validate_patient, patient)
        patient = Patient(identifier="   ")
        self.assertRaises(NeofoxDataValidationException, ModelValidator.validate_patient, patient)

    def test_bad_is_rna_available(self):
        ModelValidator.validate_patient(Patient(identifier="123", is_rna_available=True))
        ModelValidator.validate_patient(Patient(identifier="123", is_rna_available=False))
        self.assertRaises(NeofoxDataValidationException, ModelValidator.validate_patient,
                          Patient(identifier="123", is_rna_available="False"))
