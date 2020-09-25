from unittest import TestCase

from neofox.exceptions import NeofoxDataValidationException
from neofox.model.neoantigen import Gene, Neoantigen, Patient
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
