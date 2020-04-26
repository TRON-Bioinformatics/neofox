#!/usr/bin/env python

from Bio import SeqIO
from logzero import logger

import input.aa_index.aa_index as aa_index
from input import MHC_I, MHC_II
from input.epitope import Epitope
from input.helpers import data_import
from input.helpers.runner import Runner
from input.new_features import conservation_scores
from input.references import ReferenceFolder, DependenciesConfiguration


class Bunchepitopes:

    def __init__(self):
        self.references = ReferenceFolder()
        self.configuration = DependenciesConfiguration()
        self.runner = Runner()
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
        self.hlaII_available_MixMHC2pred = []
        self.rna_avail = {}

    @staticmethod
    def load_proteome(fasta_proteome):
        """
        Loads proteome in fasta format into dictionary
        """
        proteome_dict = {}
        for record in SeqIO.parse(fasta_proteome, "fasta"):
            proteome_dict[record.seq] = record.id
        return proteome_dict

    @staticmethod
    def load_rna_reference(rna_reference_file, tissue):
        """
        Loads RNA reference file and returns dictionary with mean, standard deviation and sum of expression for each gene
        """
        rna_dict = {}
        ref_tuple = data_import.import_dat_general(rna_reference_file)
        ref = data_import.get_data_from_tuple(ref_tuple)
        ref_head = data_import.get_header_from_tuple(ref_tuple)
        tissue = tissue.lower()
        head_cols = [col for col in ref_head if col.startswith(tissue)]
        cols_expr = [ref_head.index(col) for col in head_cols]
        gen_col = ref_head.index("gene")
        for ii, i in enumerate(ref):
            gene_name = i[gen_col]
            expr_values = tuple(float(i[elem]) for elem in cols_expr)
            rna_dict[gene_name] = expr_values
        return rna_dict

    @staticmethod
    def load_nmer_frequency(frequency_file):
        """
        Loads file with information of frequeny of nmers
        """
        freq_dict = {}
        with open(frequency_file) as f:
            next(f)
            for line in f:
                w = line.rstrip().split(";")
                freq_dict[w[0]] = w[1]
        return freq_dict

    def load_ucsc_ids_epitopes(self, header, list_epis):
        """
        Returns set with ucsc ids of epitopes.
        """
        set_ids = set([])
        col_ucsc = header.index("UCSC_transcript")
        col_pos = header.index("substitution")
        for ii, i in enumerate(list_epis):
            ucsc_pos = conservation_scores.add_ucsc_id_to_list(i[col_ucsc], i[col_pos])
            set_ids.add(ucsc_pos)
        return set_ids

    @staticmethod
    def load_provean_matrix(prov_matrix, epitope_ids):
        """
        Loads provean scores as dictionary, but only for ucsc ids that are in epitope list
        """
        logger.info("Starting load of PROVEAN matrix" + prov_matrix)
        provean_matrix = {}
        with open(prov_matrix) as f:
            next(f)  # skips header
            for line in f:
                w = line.rstrip().split(";")
                ucsc_id_pos = w[-1]
                if ucsc_id_pos in epitope_ids:
                    provean_matrix[ucsc_id_pos] = w
        logger.info("PROVEAN matrix loaded")
        return provean_matrix

    def load_available_hla_alleles(self, mhc=MHC_I):
        """
        loads file with available hla alllels for netmhcpan4/netmhcIIpan prediction, returns set
        :param mhc:
        :return:
        """
        if mhc == MHC_II:
            fileMHC = self.references.available_mhc_ii
        else:
            fileMHC = self.references.available_mhc_i
        set_available_mhc = set()
        with open(fileMHC) as f:
            for line in f:
                set_available_mhc.add(line.strip())
        return set_available_mhc

    def load_available_allelles_mixMHC2pred(self):
        """
        loads file with available hla alllels for MixMHC2pred prediction, returns set
        :return:
        """
        path_to_HLAII_file = self.references.alleles_list_pred
        avail_alleles = []
        with open(path_to_HLAII_file) as f:
            for line in f:
                line = line.rstrip().lstrip()
                if line:
                    if line.startswith(("L", "A")):
                        continue
                    line1 = line.split()[0]
                    if line1 is not None:
                        avail_alleles.append(line1)
        return avail_alleles

    @staticmethod
    def load_patient_hla_I_allels(path_to_hla_file):
        """
        adds hla I alleles of patients as dictionary
        :param path_to_hla_file:
        :return:
        """
        patient_alleles_dict = {}
        with open(path_to_hla_file, "r") as f:
            for line in f:
                w = line.rstrip().split(";")
                if w[0] not in patient_alleles_dict:
                    patient_alleles_dict[w[0]] = w[2:]
        logger.info("HLA-I alleles: {}".format(patient_alleles_dict))
        return patient_alleles_dict

    @staticmethod
    def load_patient_hla_II_allels(path_to_hla_file):
        """
        adds hla II alleles of patients as dictionary
        :param path_to_hla_file:
        :return:
        """
        patient_alleles_dict = {}
        with open(path_to_hla_file, "r") as f:
            for line in f:
                w = line.rstrip().split(";")
                # cheating --> overwriting hla I --> check format
                patient_alleles_dict[w[0]] = w[2:]
        logger.info("HLA-II alleles: {}".format(patient_alleles_dict))
        return patient_alleles_dict

    @staticmethod
    def load_tumor_content_dict(path_to_patient_overview):
        """
        adds tumor content of patients as dictionary
        :param path_to_patient_overview:
        :return:
        """
        tumour_content_dict = {}
        with open(path_to_patient_overview) as f:
            header = next(f)
            header = header.rstrip().split(";")
            tc_col = header.index("est. Tumor content")
            for line in f:
                w = line.rstrip().split(";")
                patient = w[0]
                patient = patient.rstrip("/")
                tumour_content = w[tc_col]
                tumour_content_dict[patient] = tumour_content

        logger.info("Tumor content: {}".format(tumour_content_dict))
        return tumour_content_dict

    @staticmethod
    def load_rna_seq_avail_dict(path_to_patient_overview):
        """
        adds info if rna seq was available as dictionary
        :param path_to_patient_overview:
        :return:
        """
        rna_avail_dict = {}
        with open(path_to_patient_overview) as f:
            header = next(f)
            header = header.rstrip().split(";")
            rna_col = header.index("rna_avail")
            for line in f:
                w = line.rstrip().split(";")
                patient = w[0]
                patient = patient.rstrip("/")
                rna_avail = w[rna_col]
                rna_avail_dict[patient] = rna_avail
        logger.info("RNA availability: {}".format(rna_avail_dict))
        return rna_avail_dict

    def write_to_file_sorted(self, d, header):
        """Transforms dictionary (property --> epitopes). To one unit (epitope) corresponding values are concentrated in one list
        and printed ';' separated."""
        features_names = []
        for key in d:
            if key not in header:
                features_names.append(key)
        features_names.sort()
        header.extend(features_names)
        print("\t".join(header))
        for i in range(len(d["mutation"])):
            z = [d[col][i] for col in header]
            print("\t".join(z))

    def initialise_properties(self, data, path_to_hla_file, tissue, tumour_content_file):
        '''adds information to Bunchepitopes class that are needed for mutated peptide sequence annotation
        '''
        self.rna_reference = self.load_rna_reference(self.references.gtex, tissue)
        freq_file1 = self.references.aa_freq_prot
        freq_file2 = self.references.four_mer_freq
        self.aa_frequency = self.load_nmer_frequency(freq_file1)
        self.fourmer_frequency = self.load_nmer_frequency(freq_file2)
        self.aa_index1_dict = aa_index.parse_aaindex1(self.references.aaindex1)
        self.aa_index2_dict = aa_index.parse_aaindex2(self.references.aaindex2)
        prov_file = self.references.prov_scores_mapped3
        self.hla_available_alleles = self.load_available_hla_alleles()
        self.hlaII_available_alleles = self.load_available_hla_alleles(mhc=MHC_II)
        self.hlaII_available_MixMHC2pred = self.load_available_allelles_mixMHC2pred()
        self.patient_hla_I_allels = self.load_patient_hla_I_allels(path_to_hla_file)
        self.patient_hla_II_allels = self.load_patient_hla_II_allels(path_to_hla_file)
        # tumour content
        if tumour_content_file != "":
            self.tumour_content = self.load_tumor_content_dict(tumour_content_file)
            self.rna_avail = self.load_rna_seq_avail_dict(tumour_content_file)
        self.ucsc_ids = self.load_ucsc_ids_epitopes(data[0], data[1])
        self.provean_matrix = self.load_provean_matrix(prov_file, self.ucsc_ids)

    def wrapper_table_add_feature_annotation(self, file, indel, path_to_hla_file, tissue, tumour_content_file):
        """ Loads epitope data (if file has been not imported to R; colnames need to be changed), adds data to class that are needed to calculate,
        calls epitope class --> determination of epitope properties,
        write to txt file
        """
        # import epitope data
        dat = data_import.import_dat_icam(file, indel)
        if "+-13_AA_(SNV)_/_-15_AA_to_STOP_(INDEL)" in dat[0]:
            dat = data_import.change_col_names(dat)
        if "mutation_found_in_proteome" not in dat[0]:
            self.proteome_dictionary = self.load_proteome(self.references.uniprot)
        # add patient id if _mut_set.txt.transcript.squish.somatic.freq is used
        if ("patient" not in dat[0]) and ("patient.id" not in dat[0]):
            try:
                patient = file.split("/")[-3]
                if "Pt" not in patient:
                    patient = file.split("/")[-1].split(".")[0]
            except IndexError:
                patient = file.split("/")[-1].split(".")[0]
            dat[0].append("patient.id")
            for ii, i in enumerate(dat[1]):
                dat[1][ii].append(str(patient))
        # initialise information needed for feature calculation
        self.initialise_properties(dat, path_to_hla_file, tissue, tumour_content_file)
        # feature calculation for each epitope
        for ii, i in enumerate(dat[1]):
            # dict for each epitope
            z = Epitope(runner=self.runner, references=self.references, configuration=self.configuration).main(
                dat[0], dat[1][ii], self.proteome_dictionary, self.rna_reference, self.aa_frequency,
                self.fourmer_frequency, self.aa_index1_dict, self.aa_index2_dict, self.provean_matrix,
                self.hla_available_alleles, self.hlaII_available_alleles, self.patient_hla_I_allels,
                self.patient_hla_II_allels, self.tumour_content, self.hlaII_available_MixMHC2pred, self.rna_avail)
            for key in z:
                if key not in self.Allepit:
                    # keys are are feautres; values: list of feature values associated with mutated peptide sequence
                    self.Allepit[key] = [z[key]]
                else:
                    self.Allepit[key].append(z[key])
        self.write_to_file_sorted(self.Allepit, dat[0])
