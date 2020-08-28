from unittest import TestCase

from neofox.annotation_resources.uniprot.uniprot import Uniprot
import pkg_resources
import neofox


class TestUniprot(TestCase):

    @classmethod
    def setUpClass(cls):
        cls.uniprot = Uniprot(pkg_resources.resource_filename(
            neofox.tests.__name__, "resources/uniprot_human_with_isoforms.first200linesfortesting.fasta"))

    def test_sequence_not_in_uniprot(self):
        result = self.uniprot.is_sequence_not_in_uniprot("NOT_IN_UNIPROT")
        self.assertEqual(True, result)

    def test_sequence_in_uniprot(self):
        result = self.uniprot.is_sequence_not_in_uniprot("LLEKVKAHEIAWLHGTI")
        self.assertEqual(False, result)


