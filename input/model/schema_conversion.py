import pandas as pd
from pandas.io.json import json_normalize
import re
import difflib
from collections import defaultdict

from input.model.neoantigen import Neoantigen, Gene, Mutation


class SchemaConverter(object):

    @staticmethod
    def validate(model):
        """
        :type model: betterproto.Message
        :return:
        """
        return model.__bytes__()

    @staticmethod
    def icam2model(icam_file):
        """
        :param icam_file: the path to an iCaM output file
        :type icam_file: str
        :rtype: list[Neoepitope]
        """
        data = pd.read_csv(icam_file, sep='\t')
        SchemaConverter._enrich_icam_table(data)
        neoantigens = []
        for _, icam_entry in data.iterrows():
            neoantigens.append(SchemaConverter._icam_entry2model(icam_entry))
        for n in neoantigens:
            SchemaConverter.validate(n)
        return neoantigens

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
    def _icam_entry2model(icam_entry):

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

