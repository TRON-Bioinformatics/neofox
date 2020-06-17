from unittest import TestCase

from input.new_features.conservation_scores import ProveanAnnotator
from input.tests.integration_tests import integration_test_tools


class TestProveanAnnotator(TestCase):

    def setUp(self):
        self.references, self.configuration = integration_test_tools.load_references()
        self.header_epitopes = ["UCSC_transcript", "substitution"]
        self.epitopes = [
            ["uc010qbo.1", "A207S"],
            ["uc001ovh.1", "A41S"],
            ["uc001ovh.1", "A40S"],
            ["uc001tzg.1", "A154S"],
            ["uc001uir.1", "A39S"],
            ["uc001yqt.1", "A701S"],
            ["uc001zrt.1", "A1520S"],
            ["uc010umw.1", "A114S"],
            ["uc010umy.1", "A7S"]]
        self.annotator = ProveanAnnotator(
            provean_file=self.references.prov_scores_mapped3, header_epitopes=self.header_epitopes,
            epitopes=self.epitopes)

    def test_provean_annotator_loading(self):
        self.assertTrue(len(self.annotator.provean_matrix) <= len(self.epitopes))
        self.assertTrue(len(self.annotator.provean_matrix) > 0)

    def test_provean_annotator(self):
        provean_annotation = self.annotator.get_provean_annotation(
            mutated_aminoacid="S", ucsc_id_position="uc001tzg_154")
        self.assertIsNotNone(provean_annotation)
        self.assertTrue(provean_annotation != "NA")
        self.assertIsNotNone(float(provean_annotation))

    def test_provean_annotator_non_existing_aminoacid(self):
        provean_annotation = self.annotator.get_provean_annotation(
            mutated_aminoacid="NO_AA", ucsc_id_position="uc001tzg_154")
        self.assertIsNotNone(provean_annotation)
        self.assertTrue(provean_annotation == "NA")

    def test_provean_annotator_non_existing_gene(self):
        provean_annotation = self.annotator.get_provean_annotation(mutated_aminoacid="S", ucsc_id_position="nope_156")
        self.assertIsNotNone(provean_annotation)
        self.assertTrue(provean_annotation == "NA")
