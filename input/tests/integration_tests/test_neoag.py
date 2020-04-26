from unittest import TestCase, SkipTest

from input.neoag.neoag_gbm_model import NeoagCalculator
from input.helpers.runner import Runner
import input.tests.integration_tests.integration_test_tools as integration_test_tools


class TestNeoantigenFitness(TestCase):

    def setUp(self):
        self.references = integration_test_tools.load_references()
        self.fastafile = integration_test_tools.create_temp_aminoacid_fasta_file()
        self.runner = Runner()

    def test_neoag(self):
        result = NeoagCalculator(runner=self.runner).wrapper_neoag(
            props={'patient': "John Doe",
                   'best_affinity_epitope_netmhcpan4': 'DDDDDDD',
                   'best_affinity_netmhcpan4': 0,
                   'best_affinity_epitope_netmhcpan4_WT': 'DDDDDDV',
                   'pos_MUT_MHCI_affinity_epi': '12345'})
        self.assertTrue(isinstance(result, str))
        self.assertTrue(float(result) > 0)

    @SkipTest
    def test_legacy(self):
        # test with ott data set
        # file = "/projects/CM01_iVAC/immunogenicity_prediction/3rd_party_solutions/MHC_prediction_netmhcpan4/testdat_ott.txt"
        # hla_file ="/projects/SUMMIT/WP1.2/Literature_Cohorts/data_analysis/cohorts/ott/icam_ott/alleles.csv"
        file = "/projects/SUMMIT/WP1.2/immunogenicity_data/ivac/input_annotation/20190328_IS_IM_withoutfeatures.txt"
        hla_file = "/projects/SUMMIT/WP1.2/immunogenicity_data/ivac/hlahd/20190916_alleles_extended.csv"
        # test inest data set
        # file = "/flash/projects/WP3/AnFranziska/AnFranziska/head_seqs.txt"
        # hla_file = "/flash/projects/WP3/AnFranziska/AnFranziska/alleles.csv"
        dat = data_import.import_dat_icam(file, False)
        if "+-13_AA_(SNV)_/_-15_AA_to_STOP_(INDEL)" in dat[0]:
            dat = data_import.change_col_names(dat)
        # available MHC alleles
        set_available_mhc = predict_all_epitopes.Bunchepitopes().load_available_hla_alleles()
        # hla allele of patients
        patient_hlaI = predict_all_epitopes.Bunchepitopes().load_patient_hla_I_allels(hla_file)
        patient_hlaII = predict_all_epitopes.Bunchepitopes().load_patient_hla_II_allels(hla_file)

        print(patient_hlaI)
        print(patient_hlaII)

        for ii, i in enumerate(dat[1]):
            if ii < 2:
                print(ii)
                dict_epi = epitope.Epitope()
                dict_epi.init_properties(dat[0], dat[1][ii])
                dict_epi.add_features(self_similarity.position_of_mutation_epitope(dict_epi.properties, MHC_I),
                                      "pos_MUT_MHCI")
                np = netmhcpan_prediction.NetmhcpanBestPrediction()
                xmer_mut = dict_epi.properties["X..13_AA_.SNV._._.15_AA_to_STOP_.INDEL."]
                tmp_fasta_file = tempfile.NamedTemporaryFile(prefix="tmp_singleseq_", suffix=".fasta", delete=False)
                tmp_fasta = tmp_fasta_file.name
                print(tmp_fasta, file=sys.stderr)
                tmp_prediction_file = tempfile.NamedTemporaryFile(prefix="netmhcpanpred_", suffix=".csv", delete=False)
                tmp_prediction = tmp_prediction_file.name
                print(tmp_prediction, file=sys.stderr)
                np.generate_fasta(dict_epi.properties, tmp_fasta, mut=True)
                alleles = np.get_hla_allels(dict_epi.properties, patient_hlaI)
                # print alleles
                np.mhc_prediction(alleles, set_available_mhc, tmp_fasta, tmp_prediction)
                dict_epi.properties["Position_Xmer_Seq"] = np.mut_position_xmer_seq(dict_epi.properties)
                preds = np.filter_binding_predictions(dict_epi.properties, tmp_prediction)
                best_epi_affinity = np.minimal_binding_score(preds, rank=False)
                dict_epi.add_features(np.add_best_epitope_info(best_epi_affinity, "Aff(nM)"),
                                      "best_affinity_netmhcpan4")
                dict_epi.add_features(np.add_best_epitope_info(best_epi_affinity, "Icore"),
                                      "best_affinity_epitope_netmhcpan4 ")
                dict_epi.add_features(np.add_best_epitope_info(best_epi_affinity, "HLA"), "best4_affinity_allele")
                xmer_wt = dict_epi.properties["X.WT._..13_AA_.SNV._._.15_AA_to_STOP_.INDEL."]
                # print >> sys.stderr, "WT seq: " + xmer_wt
                tmp_fasta_file = tempfile.NamedTemporaryFile(prefix="tmp_singleseq_", suffix=".fasta", delete=False)
                tmp_fasta = tmp_fasta_file.name
                tmp_prediction_file = tempfile.NamedTemporaryFile(prefix="netmhcpanpred_", suffix=".csv", delete=False)
                tmp_prediction = tmp_prediction_file.name
                print(tmp_prediction, file=sys.stderr)
                np = netmhcpan_prediction.NetmhcpanBestPrediction()
                np.generate_fasta(dict_epi.properties, tmp_fasta, mut=False)
                np.mhc_prediction(alleles, set_available_mhc, tmp_fasta, tmp_prediction)
                preds = np.filter_binding_predictions(dict_epi.properties, tmp_prediction)
                best_epi_affinity = np.filter_for_WT_epitope(preds,
                                                             dict_epi.properties["best_affinity_epitope_netmhcpan4"],
                                                             dict_epi.properties["best4_affinity_allele"])
                dict_epi.add_features(np.add_best_epitope_info(best_epi_affinity, "Aff(nM)"),
                                      "best_affinity_netmhcpan4_WT")
                dict_epi.add_features(np.add_best_epitope_info(best_epi_affinity, "Icore"),
                                      "best_affinity_epitope_netmhcpan4_WT")
                dict_epi.add_features(self_similarity.position_of_mutation_epitope_affinity(dict_epi.properties),
                                      "pos_MUT_MHCI_affinity_epi")

                sc = wrapper_neoag(dict_epi.properties)
                print(sc, file=sys.stderr)
                print(type(sc), file=sys.stderr)
