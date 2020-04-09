from argparse import ArgumentParser
from input.predict_all_epitopes import Bunchepitopes
import logging


def input_cli():
    parser = ArgumentParser(description='adds patient information given in sample file of a cohort to merged icam file')
    parser.add_argument('-i', '--icam_file', dest='icam_file', help='define iCaM file which should be annotated', required=True )
    parser.add_argument('-a', '--allele_file', dest='allele_file', help='define file with hla alleles of patients', required=True )
    #parser.add_argument('-r', '--reference_transcriptome', dest='ref_file', help='define suitable RNA expression reference file', default="/projects/CM27_IND_patients/GTEX_normal_tissue_data/Skin .csv" )
    parser.add_argument('-t', '--tissue', dest='tissue', help='define tissue of cancer origin', default="skin" )
    parser.add_argument('-f', '--frameshift', dest='frameshift', help='indicate by true or false if frameshift mutations or SNVs are to be considered', default=False)
    parser.add_argument('-tc', '--tumour_content', dest='tumour_content', help='pass csv file with tumour content of patient; e.g. patient_overview file ', default=False)
    args = parser.parse_args()

    icam_file = args.icam_file
    allele_file = args.allele_file
    #rna_ref = args.ref_file
    tissue = args.tissue
    indel = args.frameshift
    if args.tumour_content:
        tumour_content_file = args.tumour_content
    else:
        tumour_content_file = ""

    indel = False

    bunchepitopes = Bunchepitopes()
    logging.info("Starting INPuT...")
    bunchepitopes.wrapper_table_add_feature_annotation(icam_file, indel, allele_file, tissue, tumour_content_file)
    logging.info("Finished INPuT...")
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
