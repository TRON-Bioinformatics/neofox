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
        SchemaConverter._enrich_icam_table(data)
        neoantigens = []
        for _, icam_entry in data.iterrows():
            neoantigens.append(SchemaConverter._icam_entry2model(icam_entry, patient_id=patient_id))
        for n in neoantigens:
            SchemaConverter.validate(n)
        return neoantigens

    @staticmethod
    def patient_metadata2model(hla_file, tumor_content_file):
        """
        :param hla_file: the path to a file with the HLAs per patient for both MHC I and MHC II
        :type hla_file: str
        :param tumor_content_file: the path to a file with the tumoe content per patient
        :type tumor_content_file: str
        :rtype: list[Patient]
        """
        # parse HLA table and add metadata to patient
        alleles_stacked = SchemaConverter._parse_hlas_table(hla_file)
        patients = {}
        for patient_identifier in alleles_stacked.patient_id:
            patient = Patient()
            patient.identifier = patient_identifier
            patient.mhc_i_alleles = list(alleles_stacked[(alleles_stacked.patient_id == patient_identifier) &
                                                    (alleles_stacked.mhc_type == 'mhc_I_selection')].allele)
            patient.mhc_i_i_alleles = list(alleles_stacked[(alleles_stacked.patient_id == patient_identifier) &
                                                         (alleles_stacked.mhc_type == 'mhc_II_selection')].allele)
            patients[patient_identifier] = patient

        # parse estimated tumor content file and add metadata to patient
        tumor_content = SchemaConverter._parse_tumor_content_table(tumor_content_file)
        for patient_identifier in tumor_content.Patient:
            patient = patients.get(patient_identifier)
            if patient is not None:
                patient.estimated_tumor_content = tumor_content[tumor_content.Patient == patient_identifier][
                    'est. Tumor content'].iloc[0]
                patient.is_rna_available = tumor_content[tumor_content.Patient == patient_identifier][
                    'rna_avail'].iloc[0]

        return list(patients.values())

    @staticmethod
    def _parse_tumor_content_table(tumor_content_file):
        tumor_content = pd.read_csv(tumor_content_file, sep=';')
        tumor_content.Patient = tumor_content.Patient.transform(lambda x: x.strip('/'))
        return tumor_content

    @staticmethod
    def _parse_hlas_table(hla_file):
        alleles = pd.read_csv(hla_file, sep=';', header=None,
                              names=['patient_id', 'mhc_type'] + list(range(50))).dropna(axis=1, how='all')
        alleles_stacked = alleles.set_index(['patient_id', 'mhc_type']).stack(dropna=True).reset_index()
        del alleles_stacked['level_2']
        alleles_stacked['allele'] = alleles_stacked[0]
        del alleles_stacked[0]
        return alleles_stacked

    @staticmethod
    def model2csv(neoantigens):
        """
        :param neoantigens: list of objects of class Neoantigen
        :type neoantigens: list[Neoantigen]
        :rtype: pd.Dataframe
        """
        return json_normalize(data=[n.to_dict() for n in neoantigens])

    @staticmethod
    def csv2model(dataframe):
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
        neoantigen.clonality_estimation = None                                  # TODO: where do we get this from?
        neoantigen.expression_value = icam_entry.get('transcript_expression')   # TODO: or do we want exon expression?
        neoantigen.variant_allele_frequency = icam_entry.get('VAF_in_tumor')     # TODO: or do we want VAF in RNA?

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

