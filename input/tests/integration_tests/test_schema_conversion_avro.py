from unittest import TestCase

from input.model.avro.schema_conversion import SchemaConverter
from input.model.avro.neoantigen import Neoantigen, Gene, Mutation


class SchemaConverterTest(TestCase):

    def test_icam2model(self):
        self.icam_file = '/projects/SUMMIT/WP1.2/input/development/Pt29.sequences4testing.txt'
        # self.icam_file = '\\\\192.168.171.199\\projects$\\SUMMIT\\WP1.2\\input\\development\\Pt29.sequences4testing.txt'
        with open(self.icam_file) as f:
            self.count_lines = len(f.readlines())
        neoantigens = SchemaConverter().icam2model(self.icam_file)
        self.assertIsNotNone(neoantigens)
        self.assertEqual(self.count_lines - 1, len(neoantigens))
        for n in neoantigens:
            self.assertIsInstance(n, Neoantigen)
            g = Gene(n.gene)
            self.assertTrue(g.transcriptIdentifier is not None and len(g.transcriptIdentifier) > 0)
            m = Mutation(n.mutation)
            self.assertTrue(m.mutatedAminoacid is not None and len(m.mutatedAminoacid) > 0)