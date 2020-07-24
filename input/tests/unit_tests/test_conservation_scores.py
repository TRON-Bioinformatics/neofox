from unittest import TestCase

import pkg_resources
import input.tests
from input.new_features.conservation_scores import ProveanAnnotator


class TestProveanAnnotator(TestCase):

    def setUp(self):
        self.header_epitopes = ["UCSC_transcript", "substitution"]
        self.epitopes = [
            ["uc059atj.1", "A5S"],
            ["uc059atj.1", "A13Q"],
            ["uc059atj.1", "A27W"],
            ["uc058xwc.1", "S12Q"],
            ["uc058xwc.1", "W20S"]]
        self.annotator = ProveanAnnotator(
            provean_file=pkg_resources.resource_filename(
                input.tests.__name__, "resources/PROV_scores_mapped3.first200linesfortesting.csv"),
            header_epitopes=self.header_epitopes,
            epitopes=self.epitopes)

    def test_provean_annotator_loading(self):
        self.assertTrue(len(self.annotator.provean_matrix) <= len(self.epitopes))
        self.assertTrue(len(self.annotator.provean_matrix) > 0)

    def test_provean_annotator(self):
        provean_annotation = self.annotator.get_provean_annotation(
            mutated_aminoacid="S", ucsc_id_position="uc059atj_5")
        self.assertEqual(provean_annotation,  "-3")
        provean_annotation = self.annotator.get_provean_annotation(
            mutated_aminoacid="Q", ucsc_id_position="uc058xwc_12")
        self.assertEqual(provean_annotation,  "-4.58")

    def test_provean_annotator_non_existing_aminoacid(self):
        provean_annotation = self.annotator.get_provean_annotation(
            mutated_aminoacid="NO_AA", ucsc_id_position="uc001tzg_154")
        self.assertIsNotNone(provean_annotation)
        self.assertTrue(provean_annotation == "NA")

    def test_provean_annotator_non_existing_gene(self):
        provean_annotation = self.annotator.get_provean_annotation(mutated_aminoacid="S", ucsc_id_position="nope_156")
        self.assertIsNotNone(provean_annotation)
        self.assertTrue(provean_annotation == "NA")
