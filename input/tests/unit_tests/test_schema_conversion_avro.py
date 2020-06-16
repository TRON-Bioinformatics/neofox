import random
from unittest import TestCase

from Bio.Data import IUPACData
import json
import numpy as np

from pandas.io.json import json_normalize

from input.model.avro.schema_conversion import SchemaConverter
from input.helpers import intermediate_files
from input.model.avro.neoantigen import Neoantigen, Gene, Mutation


class SchemaConverterTest(TestCase):

    def test_model2json(self):
        neoantigens = [_get_random_neoantigen() for _ in range(5)]
        json_data = SchemaConverter.model2json(neoantigens)
        self.assertIsInstance(json_data, list)
        self.assertEqual(5, len(json_data))
        self._assert_lists_equal([Neoantigen(j) for j in json_data], neoantigens)

    def test_json2model(self):
        json_file = intermediate_files.create_temp_file(prefix='test_')
        neoantigens = [_get_random_neoantigen() for _ in range(5)]
        json_data = SchemaConverter.model2json(neoantigens)
        with open(json_file, 'w') as fd:
            fd.write(json.dumps(json_data))
        neoantigens2 = SchemaConverter.json2model(json_file)
        self.assertIsNotNone(neoantigens2)
        self._assert_lists_equal(neoantigens, neoantigens2)

    def test_model2csv(self):
        neoantigens = [_get_random_neoantigen() for _ in range(5)]
        csv_data = json_normalize(data=SchemaConverter.model2json(neoantigens))
        self.assertIsNotNone(csv_data)
        self.assertEqual(csv_data.shape[0], len(neoantigens))
        for n in neoantigens:
            self.assertEqual(n.variantAlleleFrequency,
                             csv_data[
                                 (csv_data['mutation.position'] == n.mutation['position']) &
                                 (csv_data['mutation.mutatedAminoacid'] == n.mutation['mutatedAminoacid'])
                             ].variantAlleleFrequency.iloc[0])

    def test_csv2model(self):
        neoantigens = [_get_random_neoantigen() for _ in range(5)]
        csv_data = json_normalize(data=SchemaConverter.model2json(neoantigens))
        neoantigens2 = SchemaConverter.csv2model(csv_data)
        self._assert_lists_equal(neoantigens, neoantigens2)
        
    def _assert_lists_equal(self, neoantigens, neoantigens2):
        self.assertEqual(len(neoantigens), len(neoantigens2))
        for n1, n2, in zip(neoantigens, neoantigens2):
            self.assertEqual(n1.expressionValue, n2.expressionValue)
            self.assertEqual(n1.variantAlleleFrequency, n2.variantAlleleFrequency)
            self.assertEqual(n1.clonalityEstimation, n2.clonalityEstimation)
            self.assertEqual(n1.mutation, n2.mutation)
            self.assertEqual(n1.gene, n2.gene)


class SchemaValidationTest(TestCase):
    
    def test_validation(self):
        neoantigens = [_get_random_neoantigen() for _ in range(5)]
        schema_converter = SchemaConverter()
        for n in neoantigens:
            self.assertTrue(schema_converter.validate(n))

    def test_field_not_in_the_model(self):
        neoantigen = _get_random_neoantigen()
        neoantigen._inner_dict['my_field_out_of_model'] = 5.7   # this field does not exist
        schema_converter = SchemaConverter()
        with self.assertRaises(ValueError):
            schema_converter.validate(neoantigen)

    def test_field_invalid_type(self):
        neoantigen = _get_random_neoantigen()
        neoantigen.expressionValue = "5.7"  # should be a float
        schema_converter = SchemaConverter()
        with self.assertRaises(ValueError):
            schema_converter.validate(neoantigen)

    def test_missing_non_nullable_field(self):
        neoantigen = _get_random_neoantigen()
        neoantigen.mutation = None  # non nullable field
        schema_converter = SchemaConverter()
        with self.assertRaises(ValueError):
            schema_converter.validate(neoantigen)

    def test_missing_non_nullable_nested_field(self):
        neoantigen = _get_random_neoantigen()
        mutation = Mutation()
        mutation.position = None    # non nullable field
        neoantigen.mutation = mutation
        schema_converter = SchemaConverter()
        with self.assertRaises(ValueError):
            schema_converter.validate(neoantigen)


def _get_random_neoantigen():
    neoantigen = Neoantigen()
    neoantigen.variantAlleleFrequency = np.random.uniform(0, 1)
    neoantigen.expressionValue = np.random.uniform(0, 50)
    neoantigen.clonalityEstimation = None
    mutation = Mutation()
    mutation.mutatedAminoacid = random.choices(list(IUPACData.protein_letters), k=1)[0]
    mutation.wildTypeAminoacid = random.choices(list(IUPACData.protein_letters), k=1)[0]
    mutation.position = np.random.randint(0, 1000)
    neoantigen.mutation = mutation
    gene = Gene()
    gene.gene = "BRCA2"
    gene.transcriptIdentifier = "ENST1234567"
    gene.assembly = "hg19"
    neoantigen.gene = gene
    return neoantigen
