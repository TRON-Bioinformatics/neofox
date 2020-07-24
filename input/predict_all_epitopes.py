#!/usr/bin/env python

from Bio import SeqIO
from logzero import logger

import input.aa_index.aa_index as aa_index
from input import MHC_I, MHC_II
from input.epitope import Epitope
from input.gtex.gtex import GTEx
from input.helpers import data_import
from input.helpers.properties_manager import PATIENT_ID
from input.helpers.runner import Runner
from input.new_features.conservation_scores import ProveanAnnotator
from input.references import ReferenceFolder, DependenciesConfiguration
from input.model.neoantigen import Patient


class ImmunogenicityNeoantigenPredictionToolbox:

    def __init__(self, icam_file, patient_id, patients_file):

        self.patient_id = patient_id

        self.references = ReferenceFolder()
        self.configuration = DependenciesConfiguration()
        self.runner = Runner()
        self.gtex = GTEx()
        self.hla_available_alleles = self.references.load_available_hla_alleles(mhc=MHC_I)
        self.hlaII_available_alleles = self.references.load_available_hla_alleles(mhc=MHC_II)

        # TODO: remove this empty initialisations
        self.Allepit = {}

        # import epitope data
        self.header, self.rows = data_import.import_dat_icam(icam_file)
        patients = data_import.import_patients_data(patients_file)
        if "+-13_AA_(SNV)_/_-15_AA_to_STOP_(INDEL)" in self.header:
            self.header, self.rows = data_import.change_col_names(header=self.header, data=self.rows)

        self.proteome_dictionary = self.load_proteome(self.references.uniprot)
        # adds patient to the table
        self.header.append(PATIENT_ID)
        logger.debug(self.patient_id)

        # TODO: when we are moving away from the icam table we will have a patient id for each neoantigen and
        # TODO: we will be able to pass the whole list of patients below
        self.tissue = list(filter(lambda x: x.identifier == self.patient_id, patients))[0].tissue
        for row in self.rows:
            row.append(str(self.patient_id))

        freq_file1 = self.references.aa_freq_prot
        freq_file2 = self.references.four_mer_freq
        self.aa_frequency = self.load_nmer_frequency(freq_file1)
        self.fourmer_frequency = self.load_nmer_frequency(freq_file2)
        self.aa_index1_dict = aa_index.parse_aaindex1(self.references.aaindex1)
        self.aa_index2_dict = aa_index.parse_aaindex2(self.references.aaindex2)

        self.patient_hla_I_allels = {p.identifier: p.mhc_i_alleles for p in patients}
        self.patient_hla_II_allels = {p.identifier: p.mhc_i_i_alleles for p in patients}
        self.tumour_content = {p.identifier: p.estimated_tumor_content for p in patients}
        self.rna_avail = {p.identifier: p.is_rna_available for p in patients}
        self.provean_annotator = ProveanAnnotator(provean_file=self.references.prov_scores_mapped3,
                                                  header_epitopes=self.header, epitopes=self.data)

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

    def run(self):
        """ Loads epitope data (if file has been not imported to R; colnames need to be changed), adds data to class that are needed to calculate,
        calls epitope class --> determination of epitope properties,
        write to txt file
        """
        # feature calculation for each epitope
        for row in self.rows:
            # dict for each epitope
            epitope = Epitope(
                runner=self.runner, references=self.references, configuration=self.configuration,
                provean_annotator=self.provean_annotator, gtex=self.gtex)
            features = epitope.main(
                self.header, row, self.proteome_dictionary, self.gtex, self.aa_frequency,
                self.fourmer_frequency, self.aa_index1_dict, self.aa_index2_dict,
                self.hla_available_alleles, self.hlaII_available_alleles, self.patient_hla_I_allels,
                self.patient_hla_II_allels, self.tumour_content, self.rna_avail, self.patient_id, self.tissue)
            for key in features:
                if key not in self.Allepit:
                    # keys are are feautres; values: list of feature values associated with mutated peptide sequence
                    self.Allepit[key] = [features[key]]
                else:
                    self.Allepit[key].append(features[key])
        self.write_to_file_sorted(self.Allepit, self.header)
