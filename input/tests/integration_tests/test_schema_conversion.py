from unittest import TestCase

from input.model.schema_conversion import SchemaConverter
from input.model.neoantigen import Neoantigen, Gene, Mutation, Patient


class SchemaConverterTest(TestCase):

    def test_icam2model(self):
        icam_file = '/projects/SUMMIT/WP1.2/input/development/Pt29.sequences4testing.txt'
        # self.icam_file = '\\\\192.168.171.199\\projects$\\SUMMIT\\WP1.2\\input\\development\\Pt29.sequences4testing.txt'
        with open(icam_file) as f:
            self.count_lines = len(f.readlines())
        neoantigens = SchemaConverter().icam2model(icam_file)
        self.assertIsNotNone(neoantigens)
        self.assertEqual(self.count_lines - 1, len(neoantigens))
        for n in neoantigens:
            self.assertIsInstance(n, Neoantigen)
            self.assertIsInstance(n.gene, Gene)
            self.assertIsInstance(n.mutation, Mutation)
            self.assertTrue(n.gene.transcript_identifier is not None and len(n.gene.transcript_identifier) > 0)
            self.assertTrue(n.mutation.mutated_aminoacid is not None and len(n.mutation.mutated_aminoacid) == 1)
            self.assertTrue(n.rna_variant_allele_frequency is None or
                            (0 <= n.rna_variant_allele_frequency <= 1))
            self.assertTrue(n.rna_expression is None or n.rna_expression >= 0)
            self.assertTrue(0 <= n.dna_variant_allele_frequency <= 1)

    def test_overriding_patient_id(self):
        icam_file = '/projects/SUMMIT/WP1.2/input/development/Pt29.sequences4testing.txt'
        #icam_file = '\\\\192.168.171.199\\projects$\\SUMMIT\\WP1.2\\input\\development\\Pt29.sequences4testing.txt'
        with open(icam_file) as f:
            self.count_lines = len(f.readlines())
        neoantigens = SchemaConverter().icam2model(icam_file, patient_id='patientX')
        for n in neoantigens:
            self.assertEqual(n.patient_identifier, 'patientX')
        neoantigens = SchemaConverter().icam2model(icam_file)
        for n in neoantigens:
            self.assertEqual(n.patient_identifier, None)
