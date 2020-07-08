from argparse import ArgumentParser

from logzero import logger

from input.predict_all_epitopes import BunchEpitopes


def input_cli():
    parser = ArgumentParser(description='adds patient information given in sample file of a cohort to merged icam file')
    parser.add_argument('-i', '--icam-file', dest='icam_file', help='define iCaM file which should be annotated',
                        required=True)
    # TODO: once we support the input from the models this parameter will not be required
    parser.add_argument('-p', '--patient-id', dest='patient_id', help='the patient id for the iCaM file',
                        required=True)
    parser.add_argument('-d', '--patients-data', dest='patients_data',
                        help='file with data for patients with columns: identifier, estimated_tumor_content, '
                             'is_rna_available, mhc_i_alleles, mhc_i_i_alleles, tissue',
                        required=True)
    args = parser.parse_args()

    icam_file = args.icam_file
    patient_id = args.patient_id
    patients_data = args.patients_data

    bunchepitopes = BunchEpitopes()
    logger.info("Starting INPuT...")
    bunchepitopes.wrapper_table_add_feature_annotation(icam_file, patient_id, patients_data)
    logger.info("Finished INPuT...")
    '''
    file = "/projects/CM01_iVAC/immunogenicity_prediction/3rd_party_solutions/INPuT/nonprogramm_files/test_SD.csv"
    # file = "/projects/CM01_iVAC/immunogenicity_prediction/3rd_party_solutions/INPuT/nonprogramm_files/test_fulldat.txt"
    #file = "/projects/CM01_iVAC/immunogenicity_prediction/3rd_party_solutions/20170713_IS_IM_data.complete.update_Dv10.csv.annotation.csv_v2.csv"
    # file = "/projects/SUMMIT/WP1.2/Literature_Cohorts/data_analysis/20190117_hugo_prelim_sample_annot.txt"
    #file = "/projects/SUMMIT/WP1.2/Literature_Cohorts/data_analysis/20190121_hugo_merged_dat.txt"
    indel = False
    fasta_proteome = "/projects/data/human/2018_uniprot_with_isoforms/uniprot_human_with_isoforms.fasta"
    ref_file = "/projects/CM27_IND_patients/GTEX_normal_tissue_data/Skin .csv"
    path_to_hla_file = "/projects/CM01_iVAC/immunogenicity_prediction/3rd_party_solutions/indels/RB_0004_labHLA_V2.csv"
    #path_to_hla_file = "/projects/SUMMIT/WP1.2/Literature_Cohorts/data_analysis/alleles.csv"


    Bunchepitopes().main(file, indel, fasta_proteome, ref_file, path_to_hla_file)
    '''

# def epitope_cli():
#     # file = "/projects/CM01_iVAC/immunogenicity_prediction/3rd_party_solutions/INPuT2/nonprogramm_files/test_SD.csv"
#     # file = "/projects/CM01_iVAC/immunogenicity_prediction/3rd_party_solutions/INPuT2/nonprogramm_files/test_fulldat.txt"
#     indel = False
#     fasta_proteome = "/projects/data/human/2018_uniprot_with_isoforms/uniprot_human_with_isoforms.fasta"
#     ref_file = "/projects/CM27_IND_patients/GTEX_normal_tissue_data/Skin .csv"
#     file = "/projects/CM01_iVAC/immunogenicity_prediction/3rd_party_solutions/MHC_prediction_netmhcpan4/testdat_ott.txt"
#     hla_file = "/projects/SUMMIT/WP1.2/Literature_Cohorts/data_analysis/cohorts/ott/icam_ott/alleles.csv"
#
#     # predfeatallBunchepitopes
#     predictAll = Bunchepitopes()
#     # args = parser.parse_args()
#     subprocess.call(["predict_all_epitopes", '-i', file, '-a', hla_file])
#
#     # z = Epitope().main(dat[0], dat[1][ii], self.proteome_dictionary, self.rna_reference, self.aa_frequency, self.fourmer_frequency, self.aa_index1_dict, self.aa_index2_dict, self.provean_matrix, self.hla_available_alleles, self.patient_hla_I_allels)
#
#     predictAll.main() - i
#     endTime = datetime.now()
#     print >> sys.stderr, "start: " + str(startTime) + "\nend: " + str(endTime) + "\nneeded: " + str(
#         endTime - startTime)
#     # print dat
#     # x = Epitope()
#     # x = Epitope(dat[1][1], dat[0])
#     # print vars(x)
#     # print dat[1][1][1]
#     # print dat[0][1]
#
#     # for ii,i in enumerate(dat[1]):
#     #    Epitope().main(dat[0],dat[1][ii])
#     # print x.tricks
#
#     # x.main(dat[0], dat[1][1])
#     # print x.tricks
#     # print x.tricks["transcript_position"]
#     # print dir(x)
