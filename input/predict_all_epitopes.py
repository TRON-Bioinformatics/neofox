#!/usr/bin/env python

from Bio import SeqIO
from logzero import logger

import input.aa_index.aa_index as aa_index
from input import MHC_I, MHC_II
from input.epitope import Epitope
from input.helpers import data_import
from input.helpers.properties_manager import PATIENT_ID
from input.helpers.runner import Runner
from input.new_features import conservation_scores
from input.new_features.conservation_scores import ProveanAnnotator
from input.references import ReferenceFolder, DependenciesConfiguration
from input.model.neoantigen import Patient


class BunchEpitopes:

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
        self.provean_annotator = {}
        self.hla_available_alleles = set()
        self.hlaII_available_alleles = set()
        self.patient_hla_I_alleles = {}
        self.patient_hla_II_alleles = {}
        self.tumour_content = {}
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
        ref = ref_tuple[1]
        ref_head = ref_tuple[0]
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
            z = [str(d[col][i]) for col in header]
            print("\t".join(z))

    def initialise_properties(self, header, data, patients):
        """
        adds information to Bunchepitopes class that are needed for mutated peptide sequence annotation
        :param header
        :param data:
        :type patients: list[Patient]
        :return:
        """
        # TODO: for now it reads only the tissue from the first patient in the table, we will need to do this on a per
        # TODO: patient basis
        self.rna_reference = self.load_rna_reference(self.references.gtex, patients[0].tissue)
        freq_file1 = self.references.aa_freq_prot
        freq_file2 = self.references.four_mer_freq
        self.aa_frequency = self.load_nmer_frequency(freq_file1)
        self.fourmer_frequency = self.load_nmer_frequency(freq_file2)
        self.aa_index1_dict = aa_index.parse_aaindex1(self.references.aaindex1)
        self.aa_index2_dict = aa_index.parse_aaindex2(self.references.aaindex2)
        prov_file = self.references.prov_scores_mapped3
        self.hla_available_alleles = self.references.load_available_hla_alleles(mhc=MHC_I)
        self.hlaII_available_alleles = self.references.load_available_hla_alleles(mhc=MHC_II)
        self.patient_hla_I_allels = {p.identifier: p.mhc_i_alleles for p in patients}
        self.patient_hla_II_allels = {p.identifier: p.mhc_i_i_alleles for p in patients}
        self.tumour_content = {p.identifier: p.estimated_tumor_content for p in patients}
        self.rna_avail = {p.identifier: p.is_rna_available for p in patients}
        self.provean_annotator = ProveanAnnotator(provean_file=prov_file, header_epitopes=header, epitopes=data)

    def wrapper_table_add_feature_annotation(self, icam_file, patient_id, patients_file):
        """ Loads epitope data (if file has been not imported to R; colnames need to be changed), adds data to class that are needed to calculate,
        calls epitope class --> determination of epitope properties,
        write to txt file
        """
        # import epitope data
        header, rows = data_import.import_dat_icam(icam_file)
        patients = data_import.import_patients_data(patients_file)
        if "+-13_AA_(SNV)_/_-15_AA_to_STOP_(INDEL)" in header:
            header, rows = data_import.change_col_names(header=header, data=rows)
        if "mutation_found_in_proteome" not in header:
            self.proteome_dictionary = self.load_proteome(self.references.uniprot)
        # adds patient to the table
        header.append(PATIENT_ID)
        logger.debug(patient_id)
        for row in rows:
            row.append(str(patient_id))
        # initialise information needed for feature calculation
        self.initialise_properties(header=header, data=rows, patients=patients)
        # feature calculation for each epitope
        for row in rows:
            # dict for each epitope
            epitope = Epitope(
                runner=self.runner, references=self.references, configuration=self.configuration,
                provean_annotator=self.provean_annotator)
            features = epitope.main(
                header, row, self.proteome_dictionary, self.rna_reference, self.aa_frequency,
                self.fourmer_frequency, self.aa_index1_dict, self.aa_index2_dict,
                self.hla_available_alleles, self.hlaII_available_alleles, self.patient_hla_I_allels,
                self.patient_hla_II_allels, self.tumour_content, self.rna_avail, patient_id)
            for key in features:
                if key not in self.Allepit:
                    # keys are are feautres; values: list of feature values associated with mutated peptide sequence
                    self.Allepit[key] = [features[key]]
                else:
                    self.Allepit[key].append(features[key])
        self.write_to_file_sorted(self.Allepit, header)
