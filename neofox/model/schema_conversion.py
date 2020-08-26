from itertools import islice
from typing import List

import pandas as pd
import betterproto
from betterproto import Casing
from pandas.io.json import json_normalize
import re
import difflib
from collections import defaultdict

from neofox.model.neoantigen import Neoantigen, Gene, Mutation, Patient

ICAM_FIELD_VAF_DNA = 'VAF_in_tumor'
ICAM_FIELD_VAF_RNA = 'VAF_in_RNA'
ICAM_FIELD_RNA_EXPRESSION = 'VAF_RNA_raw'
ICAM_FIELD_TRANSCRIPT_EXPRESSION = 'transcript_expression'
ICAM_FIELD_TRANSCRIPT = 'UCSC_transcript'
ICAM_FIELD_GENE = 'gene'
ICAM_FIELD_WILD_TYPE_XMER = '[WT]_+-13_AA_(SNV)_/_-15_AA_to_STOP_(INDEL)'
ICAM_FIELD_MUTATED_XMER = '+-13_AA_(SNV)_/_-15_AA_to_STOP_(INDEL)'


class SchemaConverter(object):

    @staticmethod
    def validate(model: betterproto.Message):
        # TODO: make this method capture appropriately validation issues whend dealing with int and float
        return model.__bytes__()

    @staticmethod
    def validate_neoantigen(neoantigen: Neoantigen):
        # TODO: check that all fields are consistent
        pass

    @staticmethod
    def parse_icam_file(icam_file: str, patient_id: str = None) -> List[Neoantigen]:
        """
        :param icam_file: the path to an iCaM output file
        :param patient_id: the patient identifier for all neoantigens in the iCaM file, if not provided it is
        expected as column named `patient.id` or `patient`
        :return neoantigens in model objects
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
    def parse_patients_file(patients_file: str) -> List[Patient]:
        """
        :param patients_file: the file to patients data CSV file
        :return: the parsed CSV into model objects
        """
        split_comma_separated_list = lambda x: x.split(',')
        df = pd.read_csv(
            patients_file,
            sep='\t',
            converters={'mhcIAlleles': split_comma_separated_list,
                        'mhcIIAlleles': split_comma_separated_list,
                        # TODO: remove this conversion if this is fixed
                        #  https://github.com/danielgtaylor/python-betterproto/issues/96
                        'estimatedTumorContent': lambda x: float(x) if x != "NA" else x})
        return SchemaConverter.patient_metadata_csv2objects(df)

    @staticmethod
    def objects2dataframe(model_objects: List[betterproto.Message]) -> pd.DataFrame:
        """
        :param model_objects: list of objects of subclass of betterproto.Message
        """
        return json_normalize(data=[n.to_dict(include_default_values=True) for n in model_objects])

    @staticmethod
    def object2series(model_object: betterproto.Message) -> pd.Series:
        """
        :param model_object: object of subclass of betterproto.Message
        """
        return json_normalize(data=model_object.to_dict(casing=Casing.SNAKE, include_default_values=True)).iloc[0]

    @staticmethod
    def object2flat_dict(model: betterproto.Message) -> dict:
        """
        Transforms a model object into a flat dict. Nested fields are concatenated with a dot
        """
        return SchemaConverter.object2series(model).to_dict()

    @staticmethod
    def neoantigens_csv2object(series: pd.Series) -> Neoantigen:
        """transforms an entry from a CSV into an object"""
        return Neoantigen().from_dict(SchemaConverter._flat_dict2nested_dict(flat_dict=series.to_dict()))

    @staticmethod
    def patient_metadata_csv2objects(dataframe: pd.DataFrame) -> List[Patient]:
        """transforms an patients CSV into a list of objects"""
        patients = []
        for _, row in dataframe.iterrows():
            patients.append(Patient().from_dict(row.to_dict()))
        return patients

    @staticmethod
    def _flat_dict2nested_dict(flat_dict: dict) -> dict:
        """transforms a flattened dict into a nested dict, assuming that the dot indicates a nested level"""
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
    def _icam_entry2model(icam_entry: dict, patient_id: str) -> Neoantigen:
        """parses an row from an iCaM file into a model object"""
        gene = Gene()
        gene.assembly = 'hg19'
        gene.gene = icam_entry.get(ICAM_FIELD_GENE)
        gene.transcript_identifier = icam_entry.get(ICAM_FIELD_TRANSCRIPT)

        mutation = Mutation()
        mutation.position = icam_entry.get('position')
        mutation.wild_type_xmer = icam_entry.get(ICAM_FIELD_WILD_TYPE_XMER)
        mutation.wild_type_aminoacid = icam_entry.get('wild_type_aminoacid')
        mutation.mutated_xmer = icam_entry.get(ICAM_FIELD_MUTATED_XMER)
        mutation.mutated_aminoacid = icam_entry.get('mutated_aminoacid')
        mutation.left_flanking_region = icam_entry.get('left_flanking_region')
        mutation.right_flanking_region = icam_entry.get('right_flanking_region')
        mutation.size_left_flanking_region = len(icam_entry.get('left_flanking_region'))
        mutation.size_right_flanking_region = len(icam_entry.get('right_flanking_region'))

        neoantigen = Neoantigen()
        neoantigen.patient_identifier = patient_id if patient_id else icam_entry.get('patient')
        neoantigen.mutation = mutation
        neoantigen.gene = gene
        # clonality estimation is not coming from iCaM
        neoantigen.clonality_estimation = None
        # missing RNA expression values are represented as -1
        # TODO: this is using the generic transcript expression instead of the expression read from RNA in patient
        # TODO: do we want to define some rules about when to use one or the other?
        vaf_rna_raw = icam_entry.get(ICAM_FIELD_TRANSCRIPT_EXPRESSION)
        neoantigen.rna_expression = vaf_rna_raw if vaf_rna_raw >= 0 else None
        vaf_in_rna = icam_entry.get(ICAM_FIELD_VAF_RNA)
        neoantigen.rna_variant_allele_frequency = vaf_in_rna if vaf_in_rna >= 0 else None
        neoantigen.dna_variant_allele_frequency = icam_entry.get(ICAM_FIELD_VAF_DNA)

        return neoantigen

    @staticmethod
    def _enrich_icam_table(data: pd.DataFrame):
        """parses some data from the icam table into the right fields in the model objects"""
        data['wild_type_aminoacid'] = data['substitution'].transform(lambda x: re.search("(\w)\d+\w", x).group(1))
        data['mutated_aminoacid'] = data['substitution'].transform(lambda x: re.search("\w\d+(\w)", x).group(1))
        data['position'] = data['substitution'].transform(lambda x: int(re.search("\w(\d+)\w", x).group(1)))
        data['left_flanking_region'] = data[[
            ICAM_FIELD_MUTATED_XMER, ICAM_FIELD_WILD_TYPE_XMER]].apply(
            lambda x: SchemaConverter._get_matching_region(x[0], x[1]), axis=1)
        data['right_flanking_region'] = data[[
            ICAM_FIELD_MUTATED_XMER, ICAM_FIELD_WILD_TYPE_XMER]].apply(
            lambda x: SchemaConverter._get_matching_region(x[0], x[1], match=1), axis=1)

    @staticmethod
    def _get_matching_region(sequence1: str, sequence2: str, match: int =0) -> str:
        """fetches an aligned block within two sequences"""
        match = difflib.SequenceMatcher(None, sequence1, sequence2).get_matching_blocks()[match]
        return sequence1[match.a : match.a + match.size]

