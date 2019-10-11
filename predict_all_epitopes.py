#!/usr/bin/env python

import os
import sys
from Bio import SeqIO
from argparse import ArgumentParser

from helpers import data_import
import FeatureLiterature
import epitope
from self_similarity import self_similarity
from neoantigen_fitness import neoantigen_fitness
from new_features import differential_expression
from new_features import amino_acid_frequency_scores as freq_score
from new_features import conservation_scores
from aa_index import aa_index

my_path = os.path.abspath(os.path.dirname(__file__))

class Bunchepitopes:
    def __init__(self):
        self.Allepit = {}
        self.proteome_dictionary = {}
        self.rna_reference = {}
        self.aa_frequency = {}
        self.fourmer_frequency = {}
        self.aa_index1_dict = {}
        self.aa_index2_dict = {}
        self.ucsc_ids = set([])
        self.provean_matrix = {}
        self.hla_available_alleles = set()
        self.hlaII_available_alleles = set()
        self.patient_hla_I_alleles = {}
        self.patient_hla_II_alleles = {}
        self.tumour_content = {}

    def build_proteome_dict(self, fasta_proteome):
        """Loads proteome in fasta format into dictionary
        """
        proteome_dict = {}
        for record in SeqIO.parse(fasta_proteome, "fasta"):
            proteome_dict[record.seq] = record.id
        return proteome_dict

    def add_rna_reference(self, rna_reference_file, tissue):
        """Loads RNA reference file and returns dictionary with mean, standard deviation and sum of expression for each gene
        """
        rna_dict = {}
        ref_tuple = data_import.import_dat_general(rna_reference_file)
        ref = data_import.get_data_from_tuple(ref_tuple)
        ref_head = data_import.get_header_from_tuple(ref_tuple)
        tissue = tissue.lower()
        head_cols = [col for col in ref_head if col.startswith(tissue)]
        cols_expr = [ref_head.index(col) for col in head_cols]
        gen_col = ref_head.index("gene")
        for ii,i in enumerate(ref):
            gene_name = i[gen_col]
            expr_values = tuple(float(i[elem]) for elem in cols_expr)
            rna_dict[gene_name] = expr_values
        return rna_dict

    def add_nmer_frequency(self, frequency_file):
        """Loads file with information of frequeny of nmers
        """
        freq_dict = {}
        with open(frequency_file) as f:
            header = next(f)
            for line in f:
                w = line.rstrip().split(";")
                freq_dict[w[0]] = w[1]
        return freq_dict

    def add_ucsc_ids_epitopes(self, header, list_epis):
        """Returns set with ucsc ids of epitopes.
        """
        set_ids = set([])
        col_ucsc = header.index("UCSC_transcript")
        col_pos = header.index("substitution")
        for ii,i in enumerate(list_epis):
            ucsc_pos = conservation_scores.add_ucsc_id_to_list(i[col_ucsc], i[col_pos])
            #set_ids.add(i[col_ucsc].split(".")[0])
            set_ids.add(ucsc_pos)
        return set_ids

    def add_provean_matrix(self, prov_matrix, epitope_ids):
        """Loads provean scores as dictionary, but only for ucsc ids that are in epitope list
        """
        print >> sys.stderr, "attach " + prov_matrix
        provean_matrix = {}
        with open(prov_matrix) as f:
            header = next(f)
            for line in f:
                w = line.rstrip().split(";")
                #ucsc_id = w[-2]
                ucsc_id_pos = w[-1]
                #if ucsc_id in epitope_ids:
                if ucsc_id_pos in epitope_ids:
                    provean_matrix[ucsc_id_pos] = w
        print >> sys.stderr, "attach prov matrix"
        return provean_matrix

    def add_available_hla_alleles(self, mhc="mhcI"):
        '''loads file with available hla alllels for netmhcpan4 prediction, returns set
        '''
        if mhc == "mhcII":
            fileMHC = os.path.join(my_path, "./netmhcIIpan/avail_mhcII.txt")
        else:
            fileMHC = os.path.join(my_path, "./netmhcpan4/MHC_available.csv")
        set_available_mhc = set()
        with open(fileMHC) as f:
            for line in f:
                set_available_mhc.add(line.strip())
        return set_available_mhc

    def add_patient_hla_I_allels(self, path_to_hla_file):
        '''adds hla I alleles of patients as dictionary
        '''
        patient_alleles_dict = {}
        with open(path_to_hla_file, "r") as f:
            for line in f:
                w = line.rstrip().split(";")
                if w[0] not in patient_alleles_dict:
                    patient_alleles_dict[w[0]] = w[2:]
        return patient_alleles_dict

    def add_patient_hla_II_allels(self, path_to_hla_file):
        '''adds hla II alleles of patients as dictionary
        '''
        patient_alleles_dict = {}
        with open(path_to_hla_file, "r") as f:
            for line in f:
                w = line.rstrip().split(";")
                # cheating --> overwriting hla I --> check format
                patient_alleles_dict[w[0]] = w[2:]
        return patient_alleles_dict


    def add_patient_overview(self, path_to_patient_overview):
        ''' adds tumor content of patients as dictionary
        '''
        tumour_content_dict = {}
        with open(pat_overview) as f:
            header = next(f)
            tc_col = header.index("est. Tumor content")
            for line in f:
                w = line.rstrip().split(";")
                #ucsc_id = w[-2]
                patient = w[0]
                patient = patient.rstrip("/")
                tumour_content = w[tc_col]
                tumour_content_dict[patient] = tumour_content
        return tumour_content_dict


    def write_to_file(self, d):
        """Transforms dictionary (property --> epitopes). To one unit (epitope) corresponding values are concentrated in one list
        and printed ';' separated."""
        print "\t".join(d.keys())
        for i in range(len(d["mutation"])):
            z = []
            [z.append(d[key][i]) for key in d]
            print "\t".join(z)

    def write_to_file_sorted(self, d, header):
        """Transforms dictionary (property --> epitopes). To one unit (epitope) corresponding values are concentrated in one list
        and printed ';' separated."""
        features_names = []
        for key in d:
            if key not in header:
                features_names.append(key)
        features_names.sort()
        header.extend(features_names)
        print "\t".join(header)
        for i in range(len(d["mutation"])):
            z = []
            [z.append(d[col][i]) for col in header]
            print "\t".join(z)


    def initialise_properties(self, data, db, rna_reference_file, path_to_hla_file, tissue, tumour_content_file):
        '''adds information to Bunchepitopes class that are needed for mutated peptide sequence annotation
        '''
        self.rna_reference = self.add_rna_reference(rna_reference_file, tissue)
        freq_file1 = os.path.join(my_path, "./new_features/20181108_AA_freq_prot.csv")
        freq_file2 = os.path.join(my_path, "./new_features/20181108_4mer_freq.csv")
        self.aa_frequency = self.add_nmer_frequency(freq_file1)
        self.fourmer_frequency = self.add_nmer_frequency(freq_file2)
        self.aa_index1_dict = aa_index.parse_aaindex1(os.path.join(my_path, "aa_index/aaindex1"))
        self.aa_index2_dict = aa_index.parse_aaindex2(os.path.join(my_path, "aa_index/aaindex2"))
        prov_file = os.path.join(my_path, "./new_features/PROV_scores_mapped3.csv")
        self.hla_available_alleles = self.add_available_hla_alleles()
        self.hlaII_available_alleles = self.add_available_hla_alleles(mhc = "mhcII")
        self.patient_hla_I_allels = self.add_patient_hla_I_allels(path_to_hla_file)
        self.patient_hla_II_allels = self.add_patient_hla_II_allels(path_to_hla_file)
        print >> sys.stderr, self.patient_hla_II_allels
        print >> sys.stderr, self.patient_hla_I_allels
        # tumour content
        if tumour_content_file != "":
            self.tumour_content = self.add_patient_overview(tumour_content_file)
        startTime1 = datetime.now()
        print >> sys.stderr, data[0]
        print >> sys.stderr, data[0].index("UCSC_transcript")
        self.ucsc_ids = self.add_ucsc_ids_epitopes(data[0], data[1])
        #print >> sys.stderr, self.ucsc_ids
        self.provean_matrix = self.add_provean_matrix(prov_file, self.ucsc_ids)
        endTime1 = datetime.now()
        print >> sys.stderr, "ADD PROVEAN....start: "+ str(startTime1) + "\nend: "+ str(endTime1) + "\nneeded: " + str(endTime1 - startTime1)


    def wrapper_table_add_feature_annotation(self, file, indel, db, rna_reference_file, path_to_hla_file, tissue, tumour_content_file):
        """ Loads epitope data (if file has been not imported to R; colnames need to be changed), adds data to class that are needed to calculate,
        calls epitope class --> determination of epitope properties,
        write to txt file
        """
        # import epitope data
        dat = data_import.import_dat_icam(file, indel)
        if "+-13_AA_(SNV)_/_-15_AA_to_STOP_(INDEL)" in dat[0]:
            dat = data_import.change_col_names(dat)
        if "mutation_found_in_proteome" not in dat[0]:
            self.proteome_dictionary =self.build_proteome_dict(db)
        # add patient id if _mut_set.txt.transcript.squish.somatic.freq is used
        if ("patient" not in dat[0]) and ("patient.id" not in dat[0]) :
            try:
                patient = file.split("/")[-3]
                if "Pt" not in patient:
                    patient = file.split("/")[-1].split(".")[0]
            except IndexError:
                patient = file.split("/")[-1].split(".")[0]
            dat[0].append("patient.id")
            for ii,i in enumerate(dat[1]):
                dat[1][ii].append(str(patient))
        # initialise information needed for feature calculation
        self.initialise_properties(dat, db,  rna_reference_file, path_to_hla_file, tissue, tumour_content_file)
        # feature calculation for each epitope
        for ii,i in enumerate(dat[1]):
            # dict for each epitope
            z = epitope.Epitope().main(dat[0], dat[1][ii], self.proteome_dictionary, self.rna_reference, self.aa_frequency, self.fourmer_frequency, self.aa_index1_dict, self.aa_index2_dict, self.provean_matrix, self.hla_available_alleles, self.hlaII_available_alleles, self.patient_hla_I_allels, self.patient_hla_II_allels, self.tumour_content)
            for key in z:
                if key not in self.Allepit:
                    # keys are are feautres; values: list of feature values associated with mutated peptide sequence
                    self.Allepit[key] = [z[key]]
                else:
                    self.Allepit[key].append(z[key])
        self.write_to_file_sorted(self.Allepit, dat[0])




    def main(self):
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
        db = "/projects/data/human/2018_uniprot_with_isoforms/uniprot_human_with_isoforms.fasta"
        rna_ref = "/projects/SUMMIT/WP1.2/gtex_reference/gtex_combined.csv"

        self.wrapper_table_add_feature_annotation(icam_file, indel, db, rna_ref, allele_file, tissue, tumour_content_file)





if __name__ == '__main__':
    import sys
    from datetime import datetime

    startTime = datetime.now()
    Bunchepitopes().main()




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
    endTime = datetime.now()

    print >> sys.stderr, "start: "+ str(startTime) + "\nend: "+ str(endTime) + "\nneeded: " + str(endTime - startTime)
    #print dat
    #x = Epitope()
    #x = Epitope(dat[1][1], dat[0])
    #print vars(x)
    #print dat[1][1][1]
    #print dat[0][1]


    #for ii,i in enumerate(dat[1]):
    #    Epitope().main(dat[0],dat[1][ii])
        #print x.tricks


    #x.main(dat[0], dat[1][1])
    #print x.tricks
    #print x.tricks["transcript_position"]
    #print dir(x)
