from unittest import TestCase

from input.new_features.conservation_scores import ProveanAnnotator
from input.tests.integration_tests import integration_test_tools


class TestProveanAnnotator(TestCase):

    def setUp(self):
        self.references, self.configuration = integration_test_tools.load_references()
        self.mutations = [
            "uc010qbo_207",
            "uc001ovh_41",
            "uc001ovh_40",
            "uc001tzg_154",
            "uc001uir_39",
            "uc001yqt_701",
            "uc001zrt_1520",
            "uc010umw_114",
            "uc010umy_7"]
        self.annotator = ProveanAnnotator(provean_file=self.references.prov_scores_mapped3, epitope_ids=self.mutations)

    def test_provean_annotator_loading(self):
        self.assertTrue(len(self.annotator.provean_matrix) <= len(self.mutations))
        self.assertTrue(len(self.annotator.provean_matrix) > 0)

    def test_provean_annotator(self):
        provean_annotation = self.annotator.add_provean_score_from_matrix(
            mutated_aminoacid="S", ucsc_id_position="uc001tzg_154")
        self.assertIsNotNone(provean_annotation)
        self.assertTrue(provean_annotation != "NA")
        self.assertIsNotNone(float(provean_annotation))

    def test_provean_annotator_non_existing_aminoacid(self):
        provean_annotation = self.annotator.add_provean_score_from_matrix(
            mutated_aminoacid="NO_AA", ucsc_id_position="uc001tzg_154")
        self.assertIsNotNone(provean_annotation)
        self.assertTrue(provean_annotation == "NA")

    def test_provean_annotator_non_existing_gene(self):
        provean_annotation = self.annotator.add_provean_score_from_matrix(mutated_aminoacid="S", ucsc_id_position="nope_156")
        self.assertIsNotNone(provean_annotation)
        self.assertTrue(provean_annotation == "NA")
