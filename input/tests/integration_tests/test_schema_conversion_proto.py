from unittest import TestCase

from input.model.schema_conversion import SchemaConverter
from input.model.neoantigen import Neoantigen, Gene, Mutation


class SchemaConverterTest(TestCase):

    def test_icam2model(self):
        # self.icam_file = '/projects/SUMMIT/WP1.2/input/development/Pt29.sequences4testing.txt'
        self.icam_file = '\\\\192.168.171.199\\projects$\\SUMMIT\\WP1.2\\input\\development\\Pt29.sequences4testing.txt'
        with open(self.icam_file) as f:
            self.count_lines = len(f.readlines())
        neoantigens = SchemaConverter().icam2model(self.icam_file)
        self.assertIsNotNone(neoantigens)
        self.assertEqual(self.count_lines - 1, len(neoantigens))
        for n in neoantigens:
            self.assertIsInstance(n, Neoantigen)
            self.assertIsInstance(n.gene, Gene)
            self.assertIsInstance(n.mutation, Mutation)
            self.assertTrue(n.gene.transcript_identifier is not None and len(n.gene.transcript_identifier) > 0)
            self.assertTrue(n.mutation.mutated_aminoacid is not None and len(n.mutation.mutated_aminoacid) == 1)
