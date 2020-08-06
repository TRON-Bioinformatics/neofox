from unittest import TestCase

import pkg_resources
import input.tests
from input.annotation_resources.provean.provean import ProveanAnnotator


class TestProveanAnnotator(TestCase):

    def setUp(self):
        self.annotator = ProveanAnnotator(
            provean_file=pkg_resources.resource_filename(
                input.tests.__name__, "resources/PROV_scores_mapped3.first200linesfortesting.tab.gz"))

    def test_provean_annotator_loading(self):
        self.assertIsNotNone(self.annotator.provean)
        self.assertTrue(len(self.annotator.aminoacid_indexes) > 0)

    def test_provean_annotator(self):
        provean_annotation = self.annotator.get_provean_annotation(
            mutated_aminoacid="S", protein_id="uc059atj", position=5)
        self.assertEqual(provean_annotation,  "-3")
        provean_annotation = self.annotator.get_provean_annotation(
            mutated_aminoacid="Q", protein_id="uc058xwc", position=12)
        self.assertEqual(provean_annotation,  "-4.58")

    def test_provean_annotator_non_existing_aminoacid(self):
        provean_annotation = self.annotator.get_provean_annotation(
            mutated_aminoacid="NO_AA", protein_id="uc059atj", position=5)
        self.assertIsNone(provean_annotation)

    def test_provean_annotator_non_existing_gene(self):
        provean_annotation = self.annotator.get_provean_annotation(
            mutated_aminoacid="S", protein_id="nope", position=5)
        self.assertIsNone(provean_annotation)
