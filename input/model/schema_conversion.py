from itertools import islice

import pandas as pd
from pandas.io.json import json_normalize
import re
import difflib
from collections import defaultdict

from input.model.neoantigen import Neoantigen, Gene, Mutation, Patient


class SchemaConverter(object):

    @staticmethod
    def validate(model):
        """
        :type model: betterproto.Message
        :return:
        """
        # TODO: make this method capture appropriately validation issues whend ealing with int and float
        return model.__bytes__()

    @staticmethod
    def icam2model(icam_file, patient_id=None):
        """
        :param icam_file: the path to an iCaM output file
        :type icam_file: str
        :param patient_id: the patient identifier for all neoantigens in the iCaM file, if not provided it is
        expected as column named `patient.id` or `patient`
        :type patient_id: str
        :rtype: list[Neoepitope]
        """
        data = pd.read_csv(icam_file, sep='\t')
        # filter out indels as the substitution field is reported empty by iCaM
        data = data[~data['substitution'].isna()]
        SchemaConverter._enrich_icam_table(data)
        neoantigens = []
        for _, icam_entry in data.iterrows():
            neoantigens.append(SchemaConverter._icam_entry2model(icam_entry, patient_id=patient_id))
        for n in neoantigens:
            SchemaConverter.validate(n)
        return neoantigens

    @staticmethod
    def model2csv(model_objects):
        """
        :param model_objects: list of objects of subclass of betterproto.Message
        :type model_objects: list[betterproto.Message]
        :rtype: pd.Dataframe
        """
        return json_normalize(data=[n.to_dict(include_default_values=True) for n in model_objects])

    @staticmethod
    def neoantigens_csv2model(dataframe):
        """
        :param dataframe: the input CSV in a dataframe
        :type dataframe: pd.Dataframe
        :return: the list of objects of type Neoantigen
        :rtype: list[Neoantigen]
        """
        neoantigens = []
        for _, row in dataframe.iterrows():
            neoantigens.append(Neoantigen().from_dict(SchemaConverter._flat_dict2nested_dict(flat_dict=row.to_dict())))
        return neoantigens

    @staticmethod
    def patient_metadata_csv2model(dataframe):
        """
        :param dataframe: the patient metadata CSV in a dataframe
        :type dataframe: pd.Dataframe
        :return: the list of objects of type Patient
        :rtype: list[Patient]
        """
        patients = []
        for _, row in dataframe.iterrows():
            patients.append(Patient().from_dict(row.to_dict()))
        return patients

    @staticmethod
    def _flat_dict2nested_dict(flat_dict):
        """
        :type flat_dict: dict
        :return:
        """
        nested_dict = defaultdict(lambda: {})
        for k, v in flat_dict.items():
            splitted_k = k.split('.')
            if len(splitted_k) > 2:
                raise NotImplemented("Support for dictionaries nested more than one level is not implemented")
            if len(splitted_k) == 2:
                nested_dict[splitted_k[0]][splitted_k[1]] = v
            else:
                nested_dict[k] = v
        return dict(nested_dict)

    @staticmethod
    def _icam_entry2model(icam_entry, patient_id):

        gene = Gene()
        gene.assembly = 'hg19'
        gene.gene = icam_entry.get('gene')
        gene.transcript_identifier = icam_entry.get('UCSC_transcript')

        mutation = Mutation()
        mutation.position = icam_entry.get('position')
        mutation.wild_type_aminoacid = icam_entry.get('wild_type_aminoacid')
        mutation.mutated_aminoacid = icam_entry.get('mutated_aminoacid')
        mutation.left_flanking_region = icam_entry.get('left_flanking_region')
        mutation.right_flanking_region = icam_entry.get('right_flanking_region')
        mutation.size_left_flanking_region = len(icam_entry.get('left_flanking_region'))
        mutation.size_right_flanking_region = len(icam_entry.get('right_flanking_region'))

        neoantigen = Neoantigen()
        neoantigen.patient_identifier = patient_id if patient_id else icam_entry.get('patient', icam_entry.get('patient.id'))
        neoantigen.mutation = mutation
        neoantigen.gene = gene
        # clonality estimation is not coming from iCaM
        neoantigen.clonality_estimation = None
        # missing RNA expression values are represented as -1
        vaf_rna_raw = icam_entry.get('VAF_RNA_raw')
        neoantigen.rna_expression = vaf_rna_raw if vaf_rna_raw >= 0 else None
        vaf_in_rna = icam_entry.get('VAF_in_RNA')
        neoantigen.rna_variant_allele_frequency = vaf_in_rna if vaf_in_rna >= 0 else None
        neoantigen.dna_variant_allele_frequency = icam_entry.get('VAF_in_tumor')

        return neoantigen

    @staticmethod
    def _enrich_icam_table(data):
        data['wild_type_aminoacid'] = data['substitution'].transform(lambda x: re.search("(\w)\d+\w", x).group(1))
        data['mutated_aminoacid'] = data['substitution'].transform(lambda x: re.search("\w\d+(\w)", x).group(1))
        data['position'] = data['substitution'].transform(lambda x: int(re.search("\w(\d+)\w", x).group(1)))
        data['left_flanking_region'] = data[[
            '+-13_AA_(SNV)_/_-15_AA_to_STOP_(INDEL)', '[WT]_+-13_AA_(SNV)_/_-15_AA_to_STOP_(INDEL)']].apply(
            lambda x: SchemaConverter._get_matching_region(x[0], x[1]), axis=1)
        data['right_flanking_region'] = data[[
            '+-13_AA_(SNV)_/_-15_AA_to_STOP_(INDEL)', '[WT]_+-13_AA_(SNV)_/_-15_AA_to_STOP_(INDEL)']].apply(
            lambda x: SchemaConverter._get_matching_region(x[0], x[1], match=1), axis=1)

    @staticmethod
    def _get_matching_region(sequence1, sequence2, match=0):
        match = difflib.SequenceMatcher(None, sequence1, sequence2).get_matching_blocks()[match]
        return sequence1[match.a : match.a + match.size]

