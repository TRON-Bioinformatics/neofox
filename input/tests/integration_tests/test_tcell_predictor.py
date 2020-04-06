from unittest import TestCase
from collections import defaultdict
from input.Tcell_predictor.tcellpredictor_wrapper import Tcellprediction
import input.tests.integration_tests.integration_test_tools as integration_test_tools


class TestTCellPredictor(TestCase):

    def setUp(self):
        self.references = integration_test_tools.load_references()
        self.fastafile = integration_test_tools.create_temp_aminoacid_fasta_file()

    def test_tcell_predictor(self):
        tcell_predictor = Tcellprediction(references=self.references)
        tcell_predictor.main(props=defaultdict(lambda: "blah"))
        self.assertEqual("NA", tcell_predictor.TcellPrdictionScore)
        self.assertEqual("NA", tcell_predictor.TcellPrdictionScore_9merPred)


"""
# if full icam output table is passed to script
    '''
    f = sys.argv[1]
    dat = data_import.import_dat_icam(f, indel = False)
    #print dat
    #print full_dataset(dat)
    l = full_dataset(dat, all = True)
    write_ouptut_to_file(l)
    '''

    # test for input implementation
    from input import epitope

    file = "/projects/CM01_iVAC/immunogenicity_prediction/3rd_party_solutions/MHC_prediction_netmhcpan4/testdat_ott.txt"
    dat = data_import.import_dat_icam(file, False)
    if "+-13_AA_(SNV)_/_-15_AA_to_STOP_(INDEL)" in dat[0]:
        dat = data_import.change_col_names(dat)

    path_to_Tcell_predictor = my_path

    for ii,i in enumerate(dat[1]):
        if ii < 10:
            print ii
            dict_epi = epitope.Epitope()
            dict_epi.init_properties(dat[0], dat[1][ii])
            #print dict_epi.properties
            tcellpred = Tcellprediction()

            tcellpred.main(dict_epi.properties)
            print tcellpred.TcellPrdictionScore
            print tcellpred.TcellPrdictionScore_9merPred
"""