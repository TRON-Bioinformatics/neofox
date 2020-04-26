from unittest import TestCase

import input.self_similarity.self_similarity as self_similarity
import input.tests.integration_tests.integration_test_tools as integration_test_tools
from input import MHC_I, MHC_II


class TestSelfSimilarity(TestCase):

    def setUp(self):
        self.references, self.configuration = integration_test_tools.load_references()
        self.fastafile = integration_test_tools.create_temp_aminoacid_fasta_file()

    def test_get_self_similarity_mhci(self):
        result = self_similarity.get_self_similarity(
            props={'MHC_I_epitope_.best_prediction.': 'DDD', 'MHC_I_epitope_.WT.': 'DDD'},
            mhc=MHC_I,
            references=self.references)
        self.assertEqual('1.0', result)

    def test_get_self_similarity_mhcii(self):
        result = self_similarity.get_self_similarity(
            props={'MHC_II_epitope_.best_prediction.': 'DDD', 'MHC_II_epitope_.WT.': 'DDD'},
            mhc=MHC_II,
            references=self.references)
        self.assertEqual('1.0', result)

    def test_is_improved_binder_mhci(self):
        result = self_similarity.is_improved_binder(
            props={'best%Rank_netmhcpan4': '1.0', 'best%Rank_netmhcpan4_WT': '1.3'},
            mhc=MHC_I)
        self.assertEqual('1', result)

    def test_is_improved_binder_mhcii(self):
        result = self_similarity.is_improved_binder(
            props={'MHC_II_score_.best_prediction.': '1,0', 'MHC_II_score_.WT.': '1,3'},
            mhc=MHC_II)
        self.assertEqual('1', result)
