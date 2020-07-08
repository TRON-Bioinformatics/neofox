import pandas as pd
from pandas.io.json import json_normalize
import re
import os
import difflib
import json
from collections import defaultdict

from avro_validator.schema import Schema
from input.model.avro.neoantigen import Neoantigen, Gene, Mutation

NEOANTIGEN_AVSC = 'Neoantigen.avsc'
GENE_AVSC = 'Gene.avsc'
MUTATION_AVSC = 'Mutation.avsc'


class SchemaConverter(object):

    def __init__(self):

        self.neoantigen_schema = self._initialise_schema(NEOANTIGEN_AVSC)
        self.gene_schema = self._initialise_schema(GENE_AVSC)
        self.mutation_schema = self._initialise_schema(MUTATION_AVSC)

    def _initialise_schema(self, avsc):
        schema_file = os.path.join(os.path.abspath(os.path.dirname(__file__)), avsc)
        schema = Schema(schema_file)
        return schema.parse()

    def validate(self, model):
        if isinstance(model, Neoantigen):
            valid = self.neoantigen_schema.validate(model) \
                    and self.mutation_schema.validate(model.mutation) and self.gene_schema.validate(model.gene)
        elif isinstance(model, Gene):
            valid = self.gene_schema.validate(model)
        elif isinstance(model, Mutation):
            valid = self.mutation_schema.validate(model)
        else:
            raise ValueError("Unexpected type for validation {}".format(type(model)))
        if not valid:
            raise ValueError("Invalid avro due to unknown reasons")
        return valid

    def icam2model(self, icam_file):
        """
        :param icam_file: the path to an iCaM output file
        :type icam_file: str
        :rtype: list[Neoepitope]
        """
        data = pd.read_csv(icam_file, sep='\t')
        SchemaConverter._enrich_icam_table(data)
        neoantigens = []
        for _, icam_entry in data.iterrows():
            neoantigens.append(self._icam_entry2model(icam_entry))
        for neoantigen in neoantigens:
            self.neoantigen_schema.validate(neoantigen)
        return neoantigens

    @staticmethod
    def model2json(neoantigens):
        """
        :param neoantigens: the list of objects of class Neoantigen
        :type: neoepitopes: list[Neoantigen]
        :return: the list of dict
        :rtype: list[dict]
        """
        json_data = []
        for neoantigen in neoantigens:
            temp_data = neoantigen._inner_dict
            if isinstance(neoantigen.gene, Gene):
                temp_data['gene'] = neoantigen.gene._inner_dict
            if isinstance(neoantigen.mutation, Mutation):
                temp_data['mutation'] = neoantigen.mutation._inner_dict
            json_data.append(temp_data)
        return json_data

    @staticmethod
    def json2model(json_file):
        """
        :param json_file: the path to a CSV file with the defined column names
        :type json_file: str
        :rtype: list[Neoantigen]
        """
        data = json.load(open(json_file))
        neoantigens = []
        for entry in data:
            neoantigen = Neoantigen(entry)
            neoantigen.gene = Gene(entry['gene'])
            neoantigen.mutation = Gene(entry['mutation'])
            neoantigens.append(neoantigen)
        return neoantigens

    @staticmethod
    def model2csv(neoantigens):
        """
        :param neoantigens: list of objects of class Neoantigen
        :type neoantigens: list[Neoantigen]
        :rtype: pd.Dataframe
        """
        return json_normalize(data=SchemaConverter.model2json(neoantigens))

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
            dict_data = row.to_dict()
            nested_dict = SchemaConverter._flat_dict2nested_dict(flat_dict=dict_data)
            n = Neoantigen(nested_dict)
            neoantigens.append(n)
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

    def _icam_entry2model(self, icam_entry):

        gene = Gene({
            'assembly': 'hg19',
            'gene': icam_entry.get('gene'),
            'transcriptIdentifier': icam_entry.get('UCSC_transcript')
        })

        mutation = Mutation({
            'position': icam_entry.get('position'),
            'wildTypeAminoacid': icam_entry.get('wild_type_aminoacid'),
            'mutatedAminoacid': icam_entry.get('mutated_aminoacid'),
            'leftFlankingRegion': icam_entry.get('left_flanking_region'),
            'rightFlankingRegion': icam_entry.get('right_flanking_region'),
            'sizeLeftFlankingRegion': len(icam_entry.get('left_flanking_region')),
            'sizeRightFlankingRegion': len(icam_entry.get('right_flanking_region'))
        })

        neoantigen = Neoantigen({
            'mutation': mutation._inner_dict,
            'gene': gene._inner_dict,
            'clonalityEstimation': None,                                    # TODO: where do we get this from?
            'expressionValue': icam_entry.get('transcript_expression'),     # TODO: or do we want exon expression?
            'variantAlleleFrequency': icam_entry.get('VAF_in_tumor')        # TODO: or do we want VAF in RNA?
        })
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

