#
# Copyright (c) 2020-2030 Translational Oncology at the Medical Center of the Johannes Gutenberg-University Mainz gGmbH.
#
# This file is part of Neofox
# (see https://github.com/tron-bioinformatics/neofox).
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program. If not, see <http://www.gnu.org/licenses/>.#
from typing import List

import pandas as pd
import betterproto
from betterproto import Casing
from pandas.io.json import json_normalize
import re
import difflib
from collections import defaultdict
import json
from neofox.model.neoantigen import Neoantigen, Gene, Mutation, Patient, NeoantigenAnnotations, MhcOneMolecule, \
    MhcTwoMolecule, MhcAllele, MhcOneGeneName, MhcOneGene, Zygosity, MhcTwoMoleculeName, MhcTwoGene, MhcTwoGeneName


FIELD_SUBSTITUTION = 'substitution'

FIELD_VAF_DNA = 'VAF_in_tumor'
FIELD_VAF_RNA = 'VAF_in_RNA'
FIELD_RNA_EXPRESSION = 'VAF_RNA_raw'
FIELD_TRANSCRIPT_EXPRESSION = 'transcript_expression'
FIELD_TRANSCRIPT = 'UCSC_transcript'
FIELD_GENE = 'gene'
FIELD_WILD_TYPE_XMER = '[WT]_+-13_AA_(SNV)_/_-15_AA_to_STOP_(INDEL)'
FIELD_MUTATED_XMER = '+-13_AA_(SNV)_/_-15_AA_to_STOP_(INDEL)'


class ModelConverter(object):

    HLA_ALLELE_PATTERN = re.compile(
        r"(?:HLA-)(\w+)\*?([0-9]{2}):?([0-9]{2,}):?([0-9]{2,})?:?([0-9]{2,})?([N|L|S|Q]{0,1})")

    @staticmethod
    def parse_candidate_file(candidate_file: str, patient_id: str = None) -> List[Neoantigen]:
        """
        :param candidate_file: the path to an neoantigen candidate input file
        :param patient_id: the patient identifier for all neoantigens in the input file, if not provided it is
        expected as column named `patient.id` or `patient`
        :return neoantigens in model objects
        """
        data = pd.read_csv(candidate_file, sep='\t')
        # filter out indels
        data = data[~data['substitution'].isna()]
        ModelConverter._enrich_candidate_table(data)
        neoantigens = []
        for _, candidate_entry in data.iterrows():
            neoantigens.append(ModelConverter._candidate_entry2model(candidate_entry, patient_id=patient_id))
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
                        'mhcIIAlleles': split_comma_separated_list})
        return ModelConverter.patient_metadata_csv2objects(df)

    @staticmethod
    def parse_neoantigens_file(neoantigens_file: str) -> List[Neoantigen]:
        """
        :param neoantigens_file: the file to neoantigens data CSV file
        :return: the parsed CSV into model objects
        """
        return ModelConverter.neoantigens_csv2objects(pd.read_csv(neoantigens_file, sep='\t').fillna(""))

    @staticmethod
    def objects2dataframe(model_objects: List[betterproto.Message]) -> pd.DataFrame:
        """
        :param model_objects: list of objects of subclass of betterproto.Message
        """
        return json_normalize(data=[n.to_dict(include_default_values=True) for n in model_objects])

    @staticmethod
    def objects2json(model_objects: List[betterproto.Message], output_file: str):
        """
        :param model_objects: list of objects of subclass of betterproto.Message
        """
        with open(output_file, "w") as f:
            json.dump([o.to_dict(casing=Casing.SNAKE) for o in model_objects], f)

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
        return ModelConverter.object2series(model).to_dict()

    @staticmethod
    def neoantigens_csv2object(series: pd.Series) -> Neoantigen:
        """transforms an entry from a CSV into an object"""
        return Neoantigen().from_dict(ModelConverter._flat_dict2nested_dict(flat_dict=series.to_dict()))

    @staticmethod
    def patient_metadata_csv2objects(dataframe: pd.DataFrame) -> List[Patient]:
        """transforms an patients CSV into a list of objects"""
        patients = []
        for _, row in dataframe.iterrows():
            patient_dict = row.to_dict()
            patient = Patient().from_dict(patient_dict)
            patient.mhc_one_molecules = ModelConverter._parse_mhc_one_alleles(patient_dict['mhcIAlleles'])
            patient.mhc_two_molecules = ModelConverter._parse_mhc_two_alleles(patient_dict['mhcIIAlleles'])
            patients.append(patient)
        return patients

    @staticmethod
    def neoantigens_csv2objects(dataframe: pd.DataFrame) -> List[Neoantigen]:
        """transforms an patients CSV into a list of objects"""
        neoantigens = []
        for _, row in dataframe.iterrows():
            neoantigens.append(Neoantigen().from_dict(ModelConverter._flat_dict2nested_dict(flat_dict=row.to_dict())))
        return neoantigens

    @staticmethod
    def annotations2short_wide_table(
            neoantigen_annotations: List[NeoantigenAnnotations], neoantigens: List[Neoantigen]) -> pd.DataFrame:
        dfs = []
        neoantigens_df = ModelConverter.objects2dataframe(neoantigens)
        for na in neoantigen_annotations:
            df = pd.DataFrame([a.to_dict() for a in na.annotations]).set_index('name').transpose()
            df['identifier'] = na.neoantigen_identifier
            df.reset_index(inplace=True)
            del df['index']
            dfs.append(df)
        annotations_df = pd.concat(dfs)
        return neoantigens_df.set_index('identifier').merge(annotations_df, on='identifier')

    @staticmethod
    def annotations2tall_skinny_table(neoantigen_annotations: List[NeoantigenAnnotations]) -> pd.DataFrame:
        dfs = []
        for na in neoantigen_annotations:
            df = pd.DataFrame([a.to_dict() for a in na.annotations])
            df['neoantigen_identifier'] = na.neoantigen_identifier
            dfs.append(df[df.value != "NA"])    # avoid writing NA values
        return pd.concat(dfs)

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
    def _candidate_entry2model(candidate_entry: dict, patient_id: str) -> Neoantigen:
        """parses an row from a candidate file into a model object"""
        gene = Gene()
        gene.assembly = 'hg19'
        gene.gene = candidate_entry.get(FIELD_GENE)
        gene.transcript_identifier = candidate_entry.get(FIELD_TRANSCRIPT)

        mutation = Mutation()
        mutation.position = candidate_entry.get('position')
        mutation.wild_type_aminoacid = candidate_entry.get('wild_type_aminoacid')
        mutation.mutated_aminoacid = candidate_entry.get('mutated_aminoacid')
        mutation.left_flanking_region = candidate_entry.get('left_flanking_region')
        mutation.right_flanking_region = candidate_entry.get('right_flanking_region')

        neoantigen = Neoantigen()
        neoantigen.patient_identifier = patient_id if patient_id else candidate_entry.get('patient')
        neoantigen.mutation = mutation
        neoantigen.gene = gene
        # clonality estimation is not present in candidate file at the moment
        neoantigen.clonality_estimation = None
        # missing RNA expression values are represented as -1
        vaf_rna_raw = candidate_entry.get(FIELD_TRANSCRIPT_EXPRESSION)
        neoantigen.rna_expression = vaf_rna_raw if vaf_rna_raw >= 0 else None
        vaf_in_rna = candidate_entry.get(FIELD_VAF_RNA)
        neoantigen.rna_variant_allele_frequency = vaf_in_rna if vaf_in_rna >= 0 else None
        neoantigen.dna_variant_allele_frequency = candidate_entry.get(FIELD_VAF_DNA)

        return neoantigen

    @staticmethod
    def _enrich_candidate_table(data: pd.DataFrame):
        """parses some data from the candidate table into the right fields in the model objects"""
        data['wild_type_aminoacid'] = data[FIELD_SUBSTITUTION].transform(lambda x: re.search("(\w)\d+\w", x).group(1))
        data['mutated_aminoacid'] = data[FIELD_SUBSTITUTION].transform(lambda x: re.search("\w\d+(\w)", x).group(1))
        data['position'] = data[FIELD_SUBSTITUTION].transform(lambda x: int(re.search("\w(\d+)\w", x).group(1)))
        data['left_flanking_region'] = data[[
            FIELD_MUTATED_XMER, FIELD_WILD_TYPE_XMER]].apply(
            lambda x: ModelConverter._get_matching_region(x[0], x[1]), axis=1)
        data['right_flanking_region'] = data[[
            FIELD_MUTATED_XMER, FIELD_WILD_TYPE_XMER]].apply(
            lambda x: ModelConverter._get_matching_region(x[0], x[1], match=1), axis=1)

    @staticmethod
    def _get_matching_region(sequence1: str, sequence2: str, match: int = 0) -> str:
        """fetches an aligned block within two sequences"""
        match = difflib.SequenceMatcher(None, sequence1, sequence2).get_matching_blocks()[match]
        return sequence1[match.a : match.a + match.size]

    @staticmethod
    def _parse_mhc_one_alleles(alleles: List[str]) -> List[MhcOneMolecule]:
        molecules = []
        parsed_alleles = list(map(ModelConverter.parse_mhc_allele, alleles))
        ModelConverter._validate_mhc_one_alleles(parsed_alleles)
        # do we need to validate genes anymore? add test creating MhcAllele with bad gene and see what happens
        for gene_name in MhcOneGeneName:
            gene_alleles = list(filter(lambda a: a.gene == gene_name.name, parsed_alleles))
            zygosity = ModelConverter._get_zygosity_from_alleles(gene_alleles)
            if zygosity == Zygosity.HOMOZYGOUS:
                gene_alleles = [gene_alleles[0]]   # we don't want repeated instances of the same allele
            molecules.append(MhcOneMolecule(
                name=gene_name, gene=MhcOneGene(name=gene_name, zygosity=zygosity, alleles=gene_alleles)))
        return molecules

    @staticmethod
    def _parse_mhc_two_alleles(alleles: List[str]) -> List[MhcTwoMolecule]:
        molecules = []
        parsed_alleles = list(map(ModelConverter.parse_mhc_allele, alleles))
        ModelConverter._validate_mhc_two_alleles(parsed_alleles)
        # do we need to validate genes anymore? add test creating MhcAllele with bad gene and see what happens
        for molecule_name in MhcTwoMoleculeName:
            molecule_alleles = list(filter(lambda a: molecule_name.name in a.gene, parsed_alleles))
            genes = []
            for gene_name in MhcTwoGeneName:
                gene_alleles = list(filter(lambda a: a.gene == gene_name.name, molecule_alleles))
                zygosity = ModelConverter._get_zygosity_from_alleles(gene_alleles)
                if zygosity == Zygosity.HOMOZYGOUS:
                    gene_alleles = [gene_alleles[0]]  # we don't want repeated instances of the same allele
                genes.append(MhcTwoGene(name=gene_name, zygosity=zygosity, alleles=gene_alleles))
            molecules.append(MhcTwoMolecule(name=molecule_name, genes=genes))
        return molecules

    @staticmethod
    def parse_mhc_allele(allele: str) -> MhcAllele:
        # infers gene, group and protein from the name
        match = ModelConverter.HLA_ALLELE_PATTERN.match(allele)
        assert match is not None, "Allele does not match HLA allele pattern {}".format(allele)
        gene = match.group(1)
        group = match.group(2)
        protein = match.group(3)
        name = "HLA-{gene}*{serotype}:{protein}".format(gene=gene, serotype=group, protein=protein)
        return MhcAllele(name=name, gene=gene, group=group, protein=protein)

    @staticmethod
    def _get_zygosity_from_alleles(alleles: List[MhcAllele]) -> Zygosity:
        assert len(set([a.gene for a in alleles])) <= 1, "Trying to get zygosity from alleles of different genes"
        assert len(alleles) <= 2, "More than 2 alleles for gene {}".format(alleles[0].gene)
        if len(alleles) == 2 and alleles[0].name == alleles[1].name:
            zygosity = Zygosity.HOMOZYGOUS
        elif len(alleles) == 2:
            zygosity = Zygosity.HETEROZYGOUS
        elif len(alleles) == 1:
            zygosity = Zygosity.HEMIZYGOUS
        else:
            zygosity = Zygosity.LOSS
        return zygosity

    @staticmethod
    def _validate_mhc_one_alleles(parsed_alleles: List[MhcAllele]):
        for a in parsed_alleles:
            assert a.gene in MhcOneGeneName.__members__, "Gene from MHC I allele is not valid: {}".format(a.gene)

    @staticmethod
    def _validate_mhc_two_alleles(parsed_alleles: List[MhcAllele]):
        for a in parsed_alleles:
            assert a.gene in MhcTwoGeneName.__members__, "Gene from MHC II allele is not valid: {}".format(a.gene)
