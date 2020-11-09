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
import base64
import hashlib
from typing import List, Tuple
import pandas as pd
import betterproto
import stringcase
from Bio.Alphabet.IUPAC import ExtendedIUPACProtein
from Bio.Data import IUPACData
from betterproto import Casing
from neofox.exceptions import NeofoxDataValidationException
from pandas.io.json import json_normalize
from logzero import logger
import re
import difflib
from collections import defaultdict
import json
from neofox.model.neoantigen import Neoantigen, Transcript, Mutation, Patient, NeoantigenAnnotations, Mhc2Name, \
    Mhc2GeneName, Zygosity, Mhc2Gene, Mhc2, Mhc2Isoform, MhcAllele, Mhc1Name, Mhc1, Annotation
from neofox.model.wrappers import HLA_ALLELE_PATTERN, HLA_MOLECULE_PATTERN, HLA_DR_MOLECULE_PATTERN, GENES_BY_MOLECULE, \
    get_mhc2_isoform_name
from neofox.expression_imputation import ExpressionAnnotator
from neofox.references.references import ReferenceFolder

EXTERNAL_ANNOTATIONS_NAME = "External"

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

    @staticmethod
    def parse_candidate_file(candidate_file: str, patients_dict: dict, patient_id: str = None) -> \
            Tuple[List[Neoantigen], List[NeoantigenAnnotations]]:
        """
        :param patients: a list of patient objects. this is required for imputation with gene expression
        :param candidate_file: the path to an neoantigen candidate input file
        :param patient_id: the patient identifier for all neoantigens in the input file, if not provided it is
        expected as column named `patient.id` or `patient`
        :return neoantigens in model objects + external annotations coming with the input
        """
        data = pd.read_csv(candidate_file, sep='\t')
        # filter out indels
        data = data[~data['substitution'].isna()]
        ModelConverter._enrich_candidate_table(data)
        neoantigens = []
        external_annotations = []
        for _, candidate_entry in data.iterrows():
            neoantigen = ModelConverter._candidate_entry2model(candidate_entry, patients_dict, patient_id=patient_id)
            neoantigens.append(neoantigen)
            external_annotations.append(NeoantigenAnnotations(
                neoantigen_identifier=neoantigen.identifier,
                annotations=[Annotation(name=name, value=value) for name, value in candidate_entry.iteritems()],
                annotator=EXTERNAL_ANNOTATIONS_NAME))
        return neoantigens, external_annotations

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
    def parse_neoantigens_file(neoantigens_file: str, patients_dict: dict) -> Tuple[List[Neoantigen], List[NeoantigenAnnotations]]:
        """
        :param patients_dict:
        :param neoantigens_file: the file to neoantigens data CSV file
        :return: the parsed CSV into model objects
        """
        return ModelConverter.neoantigens_csv2objects(pd.read_csv(neoantigens_file, sep='\t').fillna(""), patients_dict)

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
            patient.mhc1 = ModelConverter.parse_mhc1_alleles(patient_dict['mhcIAlleles'])
            patient.mhc2 = ModelConverter.parse_mhc2_alleles(patient_dict['mhcIIAlleles'])
            patients.append(patient)
        return patients

    @staticmethod
    def neoantigens_csv2objects(dataframe: pd.DataFrame, patients_dict: dict) -> Tuple[List[Neoantigen], List[NeoantigenAnnotations]]:
        """transforms an patients CSV into a list of objects"""
        neoantigens = []
        external_annotations = []
        for _, row in dataframe.iterrows():
            nested_dict = ModelConverter._flat_dict2nested_dict(flat_dict=row.to_dict())
            neoantigen = ModelValidator.validate_neoantigen(Neoantigen().from_dict(nested_dict))
            patient = patients_dict.get(neoantigen.patient_identifier)
            # impute gene expression if RNA-seq data is not available for a patient
            if not patient.is_rna_available:
                vaf_rna_raw = ModelConverter._substitute_expression(neoantigen, patient)
                neoantigen.rna_expression = vaf_rna_raw
            neoantigens.append(neoantigen)
            external_annotation_names = set([stringcase.snakecase(k) for k in nested_dict.keys()]).difference(
                set(Neoantigen.__annotations__.keys()))
            external_annotations.append(NeoantigenAnnotations(
                neoantigen_identifier=neoantigen.identifier,
                annotations=[Annotation(name=name, value=nested_dict.get(name)) for name in external_annotation_names],
                annotator=EXTERNAL_ANNOTATIONS_NAME))
        return neoantigens, external_annotations

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
            dfs.append(df[df.value != "NA"])  # avoid writing NA values
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
    def _candidate_entry2model(candidate_entry: dict, patients_dict: dict, patient_id: str) -> Neoantigen:
        """parses an row from a candidate file into a model object"""
        transcript = Transcript()
        transcript.assembly = 'hg19'
        transcript.gene = candidate_entry.get(FIELD_GENE)
        transcript.identifier = candidate_entry.get(FIELD_TRANSCRIPT)

        mutation = Mutation()
        mutation.position = candidate_entry.get('position')
        mutation.wild_type_aminoacid = candidate_entry.get('wild_type_aminoacid')
        mutation.mutated_aminoacid = candidate_entry.get('mutated_aminoacid')
        mutation.left_flanking_region = candidate_entry.get('left_flanking_region')
        mutation.right_flanking_region = candidate_entry.get('right_flanking_region')

        neoantigen = Neoantigen()
        neoantigen.patient_identifier = patient_id if patient_id else candidate_entry.get('patient')
        neoantigen.mutation = mutation
        neoantigen.transcript = transcript
        # clonality estimation is not present in candidate file at the moment
        neoantigen.clonality_estimation = None
        # missing RNA expression values are represented as -1
        logger.info(neoantigen.patient_identifier)
        vaf_rna_raw = candidate_entry.get(FIELD_TRANSCRIPT_EXPRESSION)
        patient = patients_dict.get(neoantigen.patient_identifier)
        # impute gene expression if RNA-seq data is not available for a patient
        if not patient.is_rna_available:
            vaf_rna_raw = ModelConverter._substitute_expression(neoantigen=neoantigen, patient=patient)
        neoantigen.rna_expression = vaf_rna_raw if vaf_rna_raw >= 0 else None
        vaf_in_rna = candidate_entry.get(FIELD_VAF_RNA)
        neoantigen.rna_variant_allele_frequency = vaf_in_rna if vaf_in_rna >= 0 else None
        neoantigen.dna_variant_allele_frequency = candidate_entry.get(FIELD_VAF_DNA)

        return ModelValidator.validate_neoantigen(neoantigen)

    @staticmethod
    def _substitute_expression(neoantigen: Neoantigen, patient: Patient) -> float:
        references = ReferenceFolder()
        expression_annotator = ExpressionAnnotator(references.tcga_expression, references.tcga_cohort_index)
        return (expression_annotator.get_gene_expression_annotation(gene_name=neoantigen.transcript.gene,
                                                                    tcga_cohort=patient.tumor_type))

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
        return sequence1[match.a: match.a + match.size]

    @staticmethod
    def parse_mhc1_alleles(alleles: List[str]) -> List[Mhc1]:
        isoforms = []
        parsed_alleles = list(map(ModelConverter.parse_mhc_allele, alleles))
        ModelConverter._validate_mhc1_alleles(parsed_alleles)
        # do we need to validate genes anymore? add test creating MhcAllele with bad gene and see what happens
        for gene_name in Mhc1Name:
            gene_alleles = list(filter(lambda a: a.gene == gene_name.name, parsed_alleles))
            zygosity = ModelConverter._get_zygosity_from_alleles(gene_alleles)
            if zygosity == Zygosity.HOMOZYGOUS:
                gene_alleles = [gene_alleles[0]]  # we don't want repeated instances of the same allele
            isoforms.append(Mhc1(name=gene_name, zygosity=zygosity, alleles=gene_alleles))
        return isoforms

    @staticmethod
    def parse_mhc2_alleles(alleles: List[str]) -> List[Mhc2]:
        mhc2s = []
        parsed_alleles = list(map(ModelConverter.parse_mhc_allele, alleles))
        ModelConverter._validate_mhc2_alleles(parsed_alleles)
        # do we need to validate genes anymore? add test creating MhcAllele with bad gene and see what happens
        for isoform_name in Mhc2Name:
            isoform_alleles = list(filter(lambda a: isoform_name.name in a.gene, parsed_alleles))
            genes = []
            for gene_name in GENES_BY_MOLECULE.get(isoform_name):
                gene_alleles = list(filter(lambda a: a.gene == gene_name.name, isoform_alleles))
                zygosity = ModelConverter._get_zygosity_from_alleles(gene_alleles)
                if zygosity == Zygosity.HOMOZYGOUS:
                    gene_alleles = [gene_alleles[0]]  # we don't want repeated instances of the same allele
                genes.append(Mhc2Gene(name=gene_name, zygosity=zygosity, alleles=gene_alleles))
            isoforms = ModelConverter._get_mhc2_isoforms(isoform_name, genes)
            mhc2s.append(Mhc2(name=isoform_name, genes=genes, isoforms=isoforms))
        return mhc2s

    @staticmethod
    def _get_mhc2_isoforms(isoform_name: Mhc2Name, genes: List[Mhc2Gene]) -> List[Mhc2Isoform]:
        isoforms = []
        if isoform_name == Mhc2Name.DR:
            assert len(genes) <= 1, "More than one gene provided for MHC II DR"
            # alpha chain of the MHC II DR is not modelled as it is constant
            isoforms = [Mhc2Isoform(name=a.name, alpha_chain=None, beta_chain=a) for g in genes for a in g.alleles]
        elif isoform_name == Mhc2Name.DP:
            assert len(genes) <= 2, "More than two genes provided for MHC II DP"
            alpha_alleles = [a for g in genes if g.name == Mhc2GeneName.DPA1 for a in g.alleles]
            beta_alleles = [a for g in genes if g.name == Mhc2GeneName.DPB1 for a in g.alleles]
            isoforms = [Mhc2Isoform(name=get_mhc2_isoform_name(a, b),
                                    alpha_chain=a, beta_chain=b) for a in alpha_alleles for b in beta_alleles]
        elif isoform_name == Mhc2Name.DQ:
            assert len(genes) <= 2, "More than two genes provided for MHC II DQ"
            alpha_alleles = [a for g in genes if g.name == Mhc2GeneName.DQA1 for a in g.alleles]
            beta_alleles = [a for g in genes if g.name == Mhc2GeneName.DQB1 for a in g.alleles]
            isoforms = [Mhc2Isoform(name=get_mhc2_isoform_name(a, b),
                                    alpha_chain=a, beta_chain=b) for a in alpha_alleles for b in beta_alleles]
        return isoforms

    @staticmethod
    def parse_mhc_allele(allele: str) -> MhcAllele:
        # infers gene, group and protein from the name
        match = HLA_ALLELE_PATTERN.match(allele)
        assert match is not None, "Allele does not match HLA allele pattern {}".format(allele)
        gene = match.group(1)
        group = match.group(2)
        protein = match.group(3)
        # builds a normalized representation of the allele
        name = "HLA-{gene}*{serotype}:{protein}".format(gene=gene, serotype=group, protein=protein)
        # ensures that full name stores the complete allele as provided but normalizes
        # its representation
        full_name = name
        six_digits_id = match.group(4)
        if six_digits_id is not None and six_digits_id != "":
            full_name = full_name + ":{}".format(six_digits_id)
            eight_digits_id = match.group(5)
            if eight_digits_id is not None and eight_digits_id != "":
                full_name = full_name + ":{}".format(eight_digits_id)
                expression_change = match.group(6)
                if expression_change is not None and expression_change != "":
                    full_name = full_name + expression_change
        return MhcAllele(full_name=full_name, name=name, gene=gene, group=group, protein=protein)

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
    def _validate_mhc1_alleles(parsed_alleles: List[MhcAllele]):
        for a in parsed_alleles:
            assert a.gene in Mhc1Name.__members__, "Gene from MHC I allele is not valid: {}".format(a.gene)

    @staticmethod
    def _validate_mhc2_alleles(parsed_alleles: List[MhcAllele]):
        for a in parsed_alleles:
            assert a.gene in Mhc2GeneName.__members__, "Gene from MHC II allele is not valid: {}".format(a.gene)


class ModelValidator(object):

    @staticmethod
    def validate(model: betterproto.Message):
        # TODO: make this method capture appropriately validation issues whend dealing with int and float
        try:
            model.__bytes__()
        except Exception as e:
            raise NeofoxDataValidationException(e)

    # TODO: add patient validation: validate GTEx tissue and MHC alleles

    @staticmethod
    def validate_neoantigen(neoantigen: Neoantigen) -> Neoantigen:

        # checks format consistency first
        ModelValidator.validate(neoantigen)

        try:
            # checks gene
            # TODO: do we want to verify existence of gene and transcript id?
            neoantigen.transcript = ModelValidator._validate_transcript(neoantigen.transcript)

            # checks mutation
            neoantigen.mutation = ModelValidator._validate_mutation(neoantigen.mutation)

            # check the expression values
            ModelValidator._validate_expression_values(neoantigen)
        except AssertionError as e:
            raise NeofoxDataValidationException(e)

        # infer other fields from the model
        return ModelValidator._enrich_neoantigen(neoantigen)

    @staticmethod
    def validate_patient(patient: Patient) -> Patient:

        # checks format consistency first
        ModelValidator.validate(patient)

        try:
            # checks that patient id is not empty considering white spaces
            patient_id = patient.identifier.strip() if patient.identifier else patient.identifier
            assert patient_id is not None and patient_id != "", "Patient identifier is empty"
            patient.identifier = patient_id

            # TODO: validate new model with isoforms, genes and alleles
            # checks MHC I
            validated_mhc1s = []
            for m in patient.mhc1:
                validated_mhc1s.append(ModelValidator._validate_mhc1(m))
            patient.mhc1 = validated_mhc1s
            # checks MHC II
            validated_mhc2s = []
            for m in patient.mhc2:
                validated_mhc2s.append(ModelValidator._validate_mhc2(m))
            patient.mhc2 = validated_mhc2s

        except AssertionError as e:
            raise NeofoxDataValidationException(e)

        return patient

    @staticmethod
    def _validate_mhc1(mhc1: Mhc1) -> Mhc1:
        assert mhc1.name in Mhc1Name, "Invalid MHC I name"
        assert mhc1.zygosity in Zygosity, "Invalid zygosity"
        alleles = mhc1.alleles
        if mhc1.zygosity in [Zygosity.HOMOZYGOUS, Zygosity.HEMIZYGOUS]:
            assert len(alleles) == 1, "An homozygous gene must have 1 allele and not {}".format(len(alleles))
        elif mhc1.zygosity == Zygosity.HETEROZYGOUS:
            assert len(alleles) == 2, "An heterozygous or hemizygous gene must have 2 alleles and not {}".format(
                len(alleles))
        elif mhc1.zygosity == Zygosity.LOSS:
            assert len(alleles) == 0, "An lost gene must have 0 alleles and not {}".format(
                len(alleles))
        validated_alleles = []
        for allele in alleles:
            validated_allele = ModelValidator.validate_mhc_allele_representation(allele)
            validated_alleles.append(validated_allele)
            assert validated_allele.gene == mhc1.name.name, \
                "The allele referring to gene {} is inside gene {}".format(validated_allele.gene, mhc1.name.name)
        mhc1.alleles = validated_alleles
        return mhc1

    @staticmethod
    def _validate_mhc2(mhc2: Mhc2) -> Mhc2:
        assert mhc2.name in Mhc2Name, "Invalid MHC II name"
        genes = mhc2.genes
        for gene in genes:
            assert gene.name in Mhc2GeneName, "Invalid gene name from MHC II"
            assert gene.name in GENES_BY_MOLECULE.get(mhc2.name), \
                "Gene {} referring to isoform {}".format(gene.name, mhc2.name.name)
            assert gene.zygosity in Zygosity, "Invalid zygosity"
            alleles = gene.alleles
            if gene.zygosity == Zygosity.HOMOZYGOUS:
                assert len(alleles) == 1, "An homozygous gene must have 1 allele and not {}".format(len(alleles))
            elif gene.zygosity in [Zygosity.HETEROZYGOUS, Zygosity.HEMIZYGOUS]:
                assert len(alleles) == 2, "An heterozygous or hemizygous gene must have 2 alleles and not {}".format(
                    len(alleles))
            elif gene.zygosity == Zygosity.LOSS:
                assert len(alleles) == 0, "An lost gene must have 0 alleles and not {}".format(
                    len(alleles))
            validated_alleles = []
            for allele in alleles:
                validated_allele = ModelValidator.validate_mhc_allele_representation(allele)
                validated_alleles.append(validated_allele)
                assert validated_allele.gene == gene.name.name, \
                    "The allele referring to gene {} is inside gene {}".format(validated_allele.gene, gene.name.name)
            gene.alleles = validated_alleles
        isoforms = mhc2.isoforms
        validated_isoforms = []
        for isoform in isoforms:
            validated_isoform = ModelValidator.validate_mhc2_isoform_representation(isoform)
            validated_isoforms.append(validated_isoform)
            if mhc2.name != Mhc2Name.DR:
                assert validated_isoform.alpha_chain.name in [a.name for g in genes for a in g.alleles], \
                    "Alpha chain allele not present in th list of alleles"
            assert validated_isoform.beta_chain.name in [a.name for g in genes for a in g.alleles], \
                "Beta chain allele not present in th list of alleles"
        mhc2.isoforms = validated_isoforms
        return mhc2

    @staticmethod
    def validate_mhc_allele_representation(allele: MhcAllele) -> MhcAllele:
        try:
            full_name = None
            if allele.full_name:
                # infers gene, group and protein from the name
                match = HLA_ALLELE_PATTERN.match(allele.full_name)
                assert match is not None, "Allele does not match HLA allele pattern {}".format(allele.name)
                gene = match.group(1)
                group = match.group(2)
                protein = match.group(3)
                full_name = allele.full_name
            elif allele.name:
                # infers gene, group and protein from the name
                match = HLA_ALLELE_PATTERN.match(allele.name)
                assert match is not None, "Allele does not match HLA allele pattern {}".format(allele.name)
                gene = match.group(1)
                group = match.group(2)
                protein = match.group(3)
            elif allele.gene and allele.group and allele.protein:
                # infers name from gene, group and protein
                gene = allele.gene
                group = allele.group
                protein = allele.protein
            else:
                raise NeofoxDataValidationException(
                    "HLA allele missing required fields, either name or gene, group and protein must be provided")

            assert gene in list(Mhc1Name.__members__.keys()) + list(Mhc2GeneName.__members__.keys()), \
                "Gene not from classic MHC: {}".format(gene)
            # builds the final allele representation and validates it just in case
            name = "HLA-{gene}*{serotype}:{protein}".format(gene=gene, serotype=group, protein=protein)
            match = HLA_ALLELE_PATTERN.match(name)
            assert match is not None, "Allele does not match HLA allele pattern {}".format(name)
        except AssertionError as e:
            raise NeofoxDataValidationException(e)

        return MhcAllele(full_name=full_name if full_name else name, name=name, gene=gene, group=group, protein=protein)

    @staticmethod
    def validate_mhc2_isoform_representation(isoform: Mhc2Isoform) -> Mhc2Isoform:
        try:
            if isoform.name:
                # infers alpha and beta chains
                match = HLA_MOLECULE_PATTERN.match(isoform.name)
                if match:
                    alpha_chain = ModelValidator.validate_mhc_allele_representation(MhcAllele(name=match.group(1)))
                    beta_chain = ModelValidator.validate_mhc_allele_representation(MhcAllele(name=match.group(2)))
                else:
                    match = HLA_DR_MOLECULE_PATTERN.match(isoform.name)
                    assert match is not None, "Molecule does not match HLA isoform pattern {}".format(isoform.name)
                    alpha_chain = MhcAllele()
                    beta_chain = ModelValidator.validate_mhc_allele_representation(MhcAllele(name=match.group(1)))
            elif isoform.alpha_chain and isoform.beta_chain:
                # infers name from gene, group and protein
                alpha_chain = ModelValidator.validate_mhc_allele_representation(isoform.alpha_chain)
                beta_chain = ModelValidator.validate_mhc_allele_representation(isoform.beta_chain)
            else:
                raise NeofoxDataValidationException("HLA isoform missing required fields")

            # builds the final allele representation and validates it just in case
            name = get_mhc2_isoform_name(alpha_chain, beta_chain)
            match = HLA_MOLECULE_PATTERN.match(name)
            match2 = HLA_DR_MOLECULE_PATTERN.match(name)
            assert match is not None or match2 is not None, \
                "Molecule does not match HLA isoform pattern {}".format(name)
        except AssertionError as e:
            raise NeofoxDataValidationException(e)

        return Mhc2Isoform(name=name, alpha_chain=alpha_chain, beta_chain=beta_chain)

    @staticmethod
    def _validate_expression_values(neoantigen):
        assert neoantigen.rna_expression is None or neoantigen.rna_expression >= 0, \
            "RNA expression should be a positive integer or zero {}".format(neoantigen.rna_expression)
        ModelValidator._validate_vaf(neoantigen.dna_variant_allele_frequency)
        ModelValidator._validate_vaf(neoantigen.rna_variant_allele_frequency)

    @staticmethod
    def _validate_mutation(mutation: Mutation) -> Mutation:
        # checks aminoacids
        mutation.mutated_aminoacid = ModelValidator._validate_aminoacid(mutation.mutated_aminoacid)
        mutation.wild_type_aminoacid = ModelValidator._validate_aminoacid(mutation.wild_type_aminoacid)

        # checks left and right flanking regions
        assert mutation.left_flanking_region is not None, "Empty left flanking region"
        mutation.left_flanking_region = mutation.left_flanking_region.strip()
        assert len(mutation.left_flanking_region) > 0, "Empty left flanking region"
        for aa in mutation.left_flanking_region:
            ModelValidator._validate_aminoacid(aa)

        assert mutation.right_flanking_region is not None, "Empty right flanking region"
        mutation.right_flanking_region = mutation.right_flanking_region.strip()
        assert len(mutation.right_flanking_region) > 0, "Empty right flanking region"
        for aa in mutation.right_flanking_region:
            ModelValidator._validate_aminoacid(aa)

        # checks the position
        assert mutation.position is not None, "Empty position"
        assert isinstance(mutation.position, int), "Position must be an integer"
        assert mutation.position > 0, "Position must be a 1-based positive integer"
        return mutation

    @staticmethod
    def _validate_transcript(transcript: Transcript) -> Transcript:

        # TODO: validate that gene symbol exists
        gene_name = transcript.gene.strip() if transcript.gene else transcript.gene
        assert gene_name is not None and len(gene_name) > 0, "Empty gene symbol"
        transcript.gene = gene_name

        # TODO: validate that transcript identifier exists
        transcript_identifier = transcript.identifier.strip() if transcript.identifier else transcript.identifier
        assert transcript_identifier is not None and len(transcript_identifier) > 0, "Empty transcript identifier"
        transcript.identifier = transcript_identifier

        # TODO: support other assemblies
        assembly = transcript.assembly if transcript.assembly else "hg19"
        assert assembly == "hg19", "Other reference genome than hg19 is not supported"
        transcript.assembly = assembly

        return transcript

    @staticmethod
    def _validate_vaf(vaf):
        assert vaf is None or 0.0 <= vaf <= 1.0, "VAF should be a positive integer or zero {}".format(vaf)

    @staticmethod
    def _enrich_neoantigen(neoantigen: Neoantigen) -> Neoantigen:
        neoantigen.mutation.wild_type_xmer = "".join([
            neoantigen.mutation.left_flanking_region,
            neoantigen.mutation.wild_type_aminoacid,
            neoantigen.mutation.right_flanking_region])
        neoantigen.mutation.mutated_xmer = "".join([
            neoantigen.mutation.left_flanking_region,
            neoantigen.mutation.mutated_aminoacid,
            neoantigen.mutation.right_flanking_region])
        neoantigen.mutation.size_left_flanking_region = len(neoantigen.mutation.left_flanking_region)
        neoantigen.mutation.size_right_flanking_region = len(neoantigen.mutation.right_flanking_region)
        neoantigen.identifier = ModelValidator.generate_neoantigen_identifier(neoantigen)
        return neoantigen

    @staticmethod
    def _validate_aminoacid(aminoacid):
        assert aminoacid is not None, "Aminoacid field cannot be empty"
        aminoacid = aminoacid.strip()
        assert isinstance(aminoacid, str), "Aminoacid has to be a string"
        if len(aminoacid) == 3:
            assert aminoacid in IUPACData.protein_letters_3to1_extended.keys(), \
                "Non existing 3 letter aminoacid {}".format(aminoacid)
            aminoacid = IUPACData.protein_letters_3to1_extended.get(aminoacid)
        if len(aminoacid) == 1:
            assert aminoacid in ExtendedIUPACProtein.letters, "Non existing aminoacid {}".format(aminoacid)
        else:
            assert False, "Invalid aminoacid {}".format(aminoacid)
        return aminoacid

    @staticmethod
    def generate_neoantigen_identifier(neoantigen: Neoantigen) -> str:
        neoantigen.identifier = None  # this needs to be done otherwise we cannot rreproduce the id after it is set
        return base64.b64encode(hashlib.md5(neoantigen.to_json().encode('utf8')).digest()).decode('utf8')
