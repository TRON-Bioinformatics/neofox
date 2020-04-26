from unittest import TestCase

import input.neoag.neoag_gbm_model as neoag_gbm_model
import input.tests.integration_tests.integration_test_tools as integration_test_tools


class TestNeoantigenFitness(TestCase):

    def setUp(self):
        self.references = integration_test_tools.load_references()
        self.fastafile = integration_test_tools.create_temp_aminoacid_fasta_file()

    def test_neoag(self):
        result = neoag_gbm_model.wrapper_neoag(
            props={'patient': "John Doe",
                   'best_affinity_epitope_netmhcpan4': 'DDDDDDD',
                   'best_affinity_netmhcpan4': 0,
                   'best_affinity_epitope_netmhcpan4_WT': 'DDDDDDV',
                   'pos_MUT_MHCI_affinity_epi': '12345'})
        self.assertTrue(isinstance(result, str))
        self.assertTrue(float(result) > 0)
