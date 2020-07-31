from unittest import TestCase

from input.annotation_resources.uniprot.uniprot import Uniprot
import pkg_resources
import input


class TestUniprot(TestCase):

    @classmethod
    def setUpClass(cls):
        cls.uniprot = Uniprot(pkg_resources.resource_filename(
            input.tests.__name__, "resources/uniprot_human_with_isoforms.first200linesfortesting.fasta"))

    def test_sequence_not_in_uniprot(self):
        result = self.uniprot.is_sequence_not_in_uniprot("NOT_IN_UNIPROT")
        self.assertEqual("1", result)

    def test_sequence_in_uniprot(self):
        result = self.uniprot.is_sequence_not_in_uniprot("LLEKVKAHEIAWLHGTI")
        self.assertEqual("0", result)


