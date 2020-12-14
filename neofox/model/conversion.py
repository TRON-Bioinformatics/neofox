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
from neofox.helpers.epitope_helper import EpitopeHelper
from neofox.exceptions import NeofoxDataValidationException
from pandas.io.json import json_normalize
from logzero import logger
from collections import defaultdict
import json
import numpy as np
from neofox.model.neoantigen import (
    Neoantigen,
    Mutation,
    Patient,
    NeoantigenAnnotations,
    Mhc2Name,
    Mhc2GeneName,
    Zygosity,
    Mhc2Gene,
    Mhc2,
    Mhc2Isoform,
    MhcAllele,
    Mhc1Name,
    Mhc1,
    Annotation,
)
from neofox.model.wrappers import (
    HLA_ALLELE_PATTERN,
    HLA_MOLECULE_PATTERN,
    HLA_DR_MOLECULE_PATTERN,
    GENES_BY_MOLECULE,
    get_mhc2_isoform_name,
)
from neofox.exceptions import NeofoxInputParametersException


EXTERNAL_ANNOTATIONS_NAME = "External"
FIELD_VAF_DNA = "VAF_in_tumor"
FIELD_VAF_RNA = "VAF_in_RNA"
FIELD_TRANSCRIPT_EXPRESSION = "transcript_expression"
FIELD_GENE = "gene"
FIELD_WILD_TYPE_XMER = "[WT]_+-13_AA_(SNV)_/_-15_AA_to_STOP_(INDEL)"
FIELD_MUTATED_XMER = "+-13_AA_(SNV)_/_-15_AA_to_STOP_(INDEL)"


class ModelConverter(object):
    @staticmethod
    def parse_candidate_file(
        candidate_file: str, patient_id: str = None
    ) -> Tuple[List[Neoantigen], List[NeoantigenAnnotations]]:
        """
        :param candidate_file: the path to an neoantigen candidate input file
        :param patient_id: the patient identifier for all neoantigens in the input file, if not provided it is
        expected as column named `patient.id` or `patient`
        :return neoantigens in model objects + external annotations coming with the input
        """
        data = pd.read_csv(candidate_file, sep="\t")
        # filter out indels
        data = data[~data["substitution"].isna()]
        data = data[~data["substitution"].str.contains("-")]
        data = data.replace({np.nan: None})
        neoantigens = []
        external_annotations = []
        for _, candidate_entry in data.iterrows():
            neoantigen = ModelConverter._candidate_entry2model(
                candidate_entry, patient_id=patient_id
            )
            neoantigens.append(neoantigen)
            external_annotations.append(
                NeoantigenAnnotations(
                    neoantigen_identifier=neoantigen.identifier,
                    annotations=[
                        # NOTE: we need to exclude the field gene from the external annotations as it matches a field
                        # in the model and thus it causes a conflict when both are renamed to gene_x and gene_y when
                        # joining
                        Annotation(name=name, value=value) for name, value
                        in candidate_entry.iteritems() if name != FIELD_GENE
                    ],
                    annotator=EXTERNAL_ANNOTATIONS_NAME,
                )
            )
        return neoantigens, external_annotations

    @staticmethod
    def parse_patients_file(patients_file: str) -> List[Patient]:
        """
        :param patients_file: the file to patients data CSV file
        :return: the parsed CSV into model objects
        """
        split_comma_separated_list = lambda x: x.split(",")
        df = pd.read_csv(
            patients_file,
            sep="\t",
            converters={
                "mhcIAlleles": split_comma_separated_list,
                "mhcIIAlleles": split_comma_separated_list,
            },
            # NOTE: forces the types of every column to avoid pandas setting the wrong type for corner cases
            dtype={
                "identifier": str
            }
        )
        return ModelConverter.patient_metadata_csv2objects(df)

    @staticmethod
    def parse_neoantigens_file(
        neoantigens_file: str,
    ) -> Tuple[List[Neoantigen], List[NeoantigenAnnotations]]:
        """
        :param neoantigens_file: the file to neoantigens data CSV file
        :return: the parsed CSV into model objects
        """
        return ModelConverter.neoantigens_csv2objects(
            pd.read_csv(
                neoantigens_file,
                sep="\t",
                # NOTE: forces the types of every column to avoid pandas setting the wrong type for corner cases
                dtype={
                    "identifier": str,
                    "gene": str,
                    "mutation.wildTypeXmer": str,
                    "mutation.mutatedXmer": str,
                    "patientIdentifier": str,
                    "dnaVariantAlleleFrequency": np.float,
                    "rnaExpression": np.float,
                    "rnaVariantAlleleFrequency": np.float
                }).replace({np.nan: None})
        )

    @staticmethod
    def parse_neoantigens_json_file(neoantigens_json_file: str) -> List[Neoantigen]:
        """
        :param neoantigens_json_file: the file to neoantigens data JSON file
        :return: the parsed JSON into model objects
        """
        return [
            Neoantigen().from_dict(n) for n in json.load(open(neoantigens_json_file))
        ]

    @staticmethod
    def objects2dataframe(model_objects: List[betterproto.Message]) -> pd.DataFrame:
        """
        :param model_objects: list of objects of subclass of betterproto.Message
        """
        return json_normalize(
            data=[n.to_dict(include_default_values=True) for n in model_objects]
        )

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
        return json_normalize(
            data=model_object.to_dict(casing=Casing.SNAKE, include_default_values=True)
        ).iloc[0]

    @staticmethod
    def object2flat_dict(model: betterproto.Message) -> dict:
        """
        Transforms a model object into a flat dict. Nested fields are concatenated with a dot
        """
        return ModelConverter.object2series(model).to_dict()

    @staticmethod
    def neoantigens_csv2object(series: pd.Series) -> Neoantigen:
        """transforms an entry from a CSV into an object"""
        return Neoantigen().from_dict(
            ModelConverter._flat_dict2nested_dict(flat_dict=series.to_dict())
        )

    @staticmethod
    def patient_metadata_csv2objects(dataframe: pd.DataFrame) -> List[Patient]:
        """transforms an patients CSV into a list of objects"""
        patients = []
        for _, row in dataframe.iterrows():
            patient_dict = row.to_dict()
            patient = Patient().from_dict(patient_dict)
            mhc_alleles = patient_dict["mhcIAlleles"]
            # NOTE: during the parsing of empty columns empty lists become a list with one empty string ...
            if len(mhc_alleles) > 1 or (len(mhc_alleles) == 1 and len(mhc_alleles[0]) > 0):
                patient.mhc1 = ModelConverter.parse_mhc1_alleles(mhc_alleles)
            else:
                patient.mhc1 = None
            mhc2_alleles = patient_dict["mhcIIAlleles"]
            if len(mhc2_alleles) > 1 or (len(mhc2_alleles) == 1 and len(mhc2_alleles[0]) > 0):
                patient.mhc2 = ModelConverter.parse_mhc2_alleles(mhc2_alleles)
            else:
                patient.mhc2 = None
            patients.append(patient)
        return patients

    @staticmethod
    def neoantigens_csv2objects(
        dataframe: pd.DataFrame,
    ) -> Tuple[List[Neoantigen], List[NeoantigenAnnotations]]:
        """transforms an patients CSV into a list of objects"""
        neoantigens = []
        external_annotations = []
        for _, row in dataframe.iterrows():
            nested_dict = ModelConverter._flat_dict2nested_dict(flat_dict=row.to_dict())
            neoantigen = ModelValidator.validate_neoantigen(Neoantigen().from_dict(nested_dict))
            neoantigens.append(neoantigen)
            external_annotation_names = set(
                [stringcase.snakecase(k) for k in nested_dict.keys()]
            ).difference(set(Neoantigen.__annotations__.keys()))
            external_annotations.append(
                NeoantigenAnnotations(
                    neoantigen_identifier=neoantigen.identifier,
                    annotations=[
                        Annotation(name=name, value=nested_dict.get(name))
                        for name in external_annotation_names
                    ],
                    annotator=EXTERNAL_ANNOTATIONS_NAME,
                )
            )
        return neoantigens, external_annotations

    @staticmethod
    def annotations2short_wide_table(
        neoantigen_annotations: List[NeoantigenAnnotations],
        neoantigens: List[Neoantigen],
    ) -> pd.DataFrame:
        dfs = []
        neoantigens_df = ModelConverter.neoantigens2table(neoantigens)
        for na in neoantigen_annotations:
            df = (
                pd.DataFrame([a.to_dict() for a in na.annotations])
                .set_index("name")
                .transpose()
            )
            df["identifier"] = na.neoantigen_identifier
            df.reset_index(inplace=True)
            del df["index"]
            dfs.append(df)
        annotations_df = pd.concat(dfs, sort=True)
        return neoantigens_df.set_index("identifier").merge(
            annotations_df, on="identifier"
        )

    @staticmethod
    def neoantigens2table(neoantigens: List[Neoantigen]) -> pd.DataFrame:
        df = ModelConverter.objects2dataframe(neoantigens)
        df["mutation.position"] = df["mutation.position"].transform(
            lambda x: ",".join([str(y) for y in x]) if x is not None else x)
        return df

    @staticmethod
    def annotations2tall_skinny_table(
        neoantigen_annotations: List[NeoantigenAnnotations],
    ) -> pd.DataFrame:
        dfs = []
        for na in neoantigen_annotations:
            df = pd.DataFrame([a.to_dict() for a in na.annotations])
            df["neoantigen_identifier"] = na.neoantigen_identifier
            dfs.append(df)  # avoid writing NA values
        return pd.concat(dfs, sort=True)

    @staticmethod
    def _flat_dict2nested_dict(flat_dict: dict) -> dict:
        """transforms a flattened dict into a nested dict, assuming that the dot indicates a nested level"""
        nested_dict = defaultdict(lambda: {})
        for k, v in flat_dict.items():
            splitted_k = k.split(".")
            if len(splitted_k) > 2:
                raise NotImplementedError(
                    "Support for dictionaries nested more than one level is not implemented"
                )
            if len(splitted_k) == 2:
                nested_dict[splitted_k[0]][splitted_k[1]] = v
            else:
                nested_dict[k] = v
        return dict(nested_dict)

    @staticmethod
    def _candidate_entry2model(candidate_entry: dict, patient_id: str) -> Neoantigen:
        """parses an row from a candidate file into a model object"""

        mutation = Mutation()
        mutation.wild_type_xmer = candidate_entry.get(FIELD_WILD_TYPE_XMER)
        mutation.mutated_xmer = candidate_entry.get(FIELD_MUTATED_XMER)

        neoantigen = Neoantigen()
        neoantigen.patient_identifier = (
            patient_id if patient_id else candidate_entry.get("patient")
        )
        if neoantigen.patient_identifier is None:
            raise NeofoxInputParametersException(
                "Please, define the parameter `patient_id` or provide a column ´patient´ in the candidate file "
            )
        neoantigen.mutation = mutation
        neoantigen.gene = candidate_entry.get(FIELD_GENE)
        # missing RNA expression values are represented as -1
        logger.info(neoantigen.patient_identifier)
        vaf_rna_raw = candidate_entry.get(FIELD_TRANSCRIPT_EXPRESSION)
        neoantigen.rna_expression = vaf_rna_raw if vaf_rna_raw >= 0 else None
        neoantigen.rna_variant_allele_frequency = candidate_entry.get(FIELD_VAF_RNA)
        neoantigen.dna_variant_allele_frequency = candidate_entry.get(FIELD_VAF_DNA)

        return ModelValidator.validate_neoantigen(neoantigen)

    @staticmethod
    def parse_mhc1_alleles(alleles: List[str]) -> List[Mhc1]:
        isoforms = []
        try:
            parsed_alleles = list(map(ModelConverter.parse_mhc_allele, alleles))
            ModelConverter._validate_mhc1_alleles(parsed_alleles)
            # do we need to validate genes anymore? add test creating MhcAllele with bad gene and see what happens
            for gene_name in Mhc1Name:
                gene_alleles = list(
                    filter(lambda a: a.gene == gene_name.name, parsed_alleles)
                )
                zygosity = ModelConverter._get_zygosity_from_alleles(gene_alleles)
                if zygosity == Zygosity.HOMOZYGOUS:
                    gene_alleles = [
                        gene_alleles[0]
                    ]  # we don't want repeated instances of the same allele
                isoforms.append(
                    Mhc1(name=gene_name, zygosity=zygosity, alleles=gene_alleles)
                )
        except AssertionError as e:
            raise NeofoxDataValidationException(e)
        return isoforms

    @staticmethod
    def parse_mhc2_alleles(alleles: List[str]) -> List[Mhc2]:
        mhc2s = []
        try:
            parsed_alleles = list(map(ModelConverter.parse_mhc_allele, alleles))
            ModelConverter._validate_mhc2_alleles(parsed_alleles)
            # do we need to validate genes anymore? add test creating MhcAllele with bad gene and see what happens
            for isoform_name in Mhc2Name:
                isoform_alleles = list(
                    filter(lambda a: isoform_name.name in a.gene, parsed_alleles)
                )
                genes = []
                for gene_name in GENES_BY_MOLECULE.get(isoform_name):
                    gene_alleles = list(
                        filter(lambda a: a.gene == gene_name.name, isoform_alleles)
                    )
                    zygosity = ModelConverter._get_zygosity_from_alleles(gene_alleles)
                    if zygosity == Zygosity.HOMOZYGOUS:
                        gene_alleles = [
                            gene_alleles[0]
                        ]  # we don't want repeated instances of the same allele
                    genes.append(
                        Mhc2Gene(
                            name=gene_name, zygosity=zygosity, alleles=gene_alleles
                        )
                    )
                isoforms = ModelConverter._get_mhc2_isoforms(isoform_name, genes)
                mhc2s.append(Mhc2(name=isoform_name, genes=genes, isoforms=isoforms))
        except AssertionError as e:
            raise NeofoxDataValidationException(e)
        return mhc2s

    @staticmethod
    def _get_mhc2_isoforms(
        isoform_name: Mhc2Name, genes: List[Mhc2Gene]
    ) -> List[Mhc2Isoform]:
        isoforms = []
        if isoform_name == Mhc2Name.DR:
            assert len(genes) <= 1, "More than one gene provided for MHC II DR"
            # alpha chain of the MHC II DR is not modelled as it is constant
            isoforms = [
                Mhc2Isoform(name=a.name, alpha_chain=MhcAllele(), beta_chain=a)
                for g in genes
                for a in g.alleles
            ]
        elif isoform_name == Mhc2Name.DP:
            assert len(genes) <= 2, "More than two genes provided for MHC II DP"
            alpha_alleles = [
                a for g in genes if g.name == Mhc2GeneName.DPA1 for a in g.alleles
            ]
            beta_alleles = [
                a for g in genes if g.name == Mhc2GeneName.DPB1 for a in g.alleles
            ]
            isoforms = [
                Mhc2Isoform(
                    name=get_mhc2_isoform_name(a, b), alpha_chain=a, beta_chain=b
                )
                for a in alpha_alleles
                for b in beta_alleles
            ]
        elif isoform_name == Mhc2Name.DQ:
            assert len(genes) <= 2, "More than two genes provided for MHC II DQ"
            alpha_alleles = [
                a for g in genes if g.name == Mhc2GeneName.DQA1 for a in g.alleles
            ]
            beta_alleles = [
                a for g in genes if g.name == Mhc2GeneName.DQB1 for a in g.alleles
            ]
            isoforms = [
                Mhc2Isoform(
                    name=get_mhc2_isoform_name(a, b), alpha_chain=a, beta_chain=b
                )
                for a in alpha_alleles
                for b in beta_alleles
            ]
        return isoforms

    @staticmethod
    def parse_mhc_allele(allele: str) -> MhcAllele:
        # infers gene, group and protein from the name
        match = HLA_ALLELE_PATTERN.match(allele)
        assert match is not None, "Allele does not match HLA allele pattern {}".format(
            allele
        )
        gene = match.group(1)
        group = match.group(2)
        protein = match.group(3)
        # builds a normalized representation of the allele
        name = "HLA-{gene}*{serotype}:{protein}".format(
            gene=gene, serotype=group, protein=protein
        )
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
        return MhcAllele(
            full_name=full_name, name=name, gene=gene, group=group, protein=protein
        )

    @staticmethod
    def parse_mhc2_isoform(isoform: str) -> Mhc2Isoform:
        # infers gene, group and protein from the name
        match = HLA_MOLECULE_PATTERN.match(isoform)
        if match:
            alpha_chain = ModelValidator.validate_mhc_allele_representation(
                MhcAllele(name=match.group(1))
            )
            beta_chain = ModelValidator.validate_mhc_allele_representation(
                MhcAllele(name=match.group(2))
            )
        else:
            match = HLA_DR_MOLECULE_PATTERN.match(isoform)
            assert (
                match is not None
            ), "Molecule does not match HLA isoform pattern {}".format(isoform)
            alpha_chain = MhcAllele()
            beta_chain = ModelValidator.validate_mhc_allele_representation(
                MhcAllele(name=match.group(1))
            )
        # builds the final allele representation and validates it just in case
        name = get_mhc2_isoform_name(alpha_chain, beta_chain)
        return Mhc2Isoform(name=name, alpha_chain=alpha_chain, beta_chain=beta_chain)

    @staticmethod
    def _get_zygosity_from_alleles(alleles: List[MhcAllele]) -> Zygosity:
        assert (
            len(set([a.gene for a in alleles])) <= 1
        ), "Trying to get zygosity from alleles of different genes"
        assert len(alleles) <= 2, "More than 2 alleles for gene {}".format(
            alleles[0].gene
        )
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
            assert (
                a.gene in Mhc1Name.__members__
            ), "Gene from MHC I allele is not valid {} at {}".format(a.gene, a.full_name)

    @staticmethod
    def _validate_mhc2_alleles(parsed_alleles: List[MhcAllele]):
        for a in parsed_alleles:
            assert (
                a.gene in Mhc2GeneName.__members__
            ), "Gene from MHC II allele is not valid {} at {}".format(a.gene, a.full_name)


class ModelValidator(object):
    @staticmethod
    def validate(model: betterproto.Message):
        # TODO: make this method capture appropriately validation issues when dealing with int and float
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
            assert neoantigen.patient_identifier is not None and len(neoantigen.patient_identifier) > 0, \
                "Missing patient identifier on neoantigen"

            # checks mutation
            neoantigen.mutation = ModelValidator._validate_mutation(neoantigen.mutation)

            # check the expression values
            ModelValidator._validate_expression_values(neoantigen)
        except AssertionError as e:
            raise NeofoxDataValidationException(e)

        # calculates the identifier now once the object is valid
        neoantigen.identifier = ModelValidator.generate_neoantigen_identifier(
            neoantigen
        )

        return neoantigen

    @staticmethod
    def validate_patient(patient: Patient) -> Patient:

        # checks format consistency first
        ModelValidator.validate(patient)

        try:
            # checks that patient id is not empty considering white spaces
            patient_id = (
                patient.identifier.strip() if patient.identifier else patient.identifier
            )
            assert (
                patient_id is not None and patient_id != ""
            ), "Patient identifier is empty"
            patient.identifier = patient_id

            # TODO: validate new model with isoforms, genes and alleles
            # checks MHC I
            validated_mhc1s = []
            if patient.mhc1:
                for m in patient.mhc1:
                    validated_mhc1s.append(ModelValidator._validate_mhc1(m))
            patient.mhc1 = validated_mhc1s
            # checks MHC II
            validated_mhc2s = []
            if patient.mhc2:
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
            assert (
                len(alleles) == 1
            ), "An homozygous gene must have 1 allele and not {}".format(len(alleles))
        elif mhc1.zygosity == Zygosity.HETEROZYGOUS:
            assert (
                len(alleles) == 2
            ), "An heterozygous or hemizygous gene must have 2 alleles and not {}".format(
                len(alleles)
            )
        elif mhc1.zygosity == Zygosity.LOSS:
            assert (
                len(alleles) == 0
            ), "An lost gene must have 0 alleles and not {}".format(len(alleles))
        validated_alleles = []
        for allele in alleles:
            validated_allele = ModelValidator.validate_mhc_allele_representation(allele)
            validated_alleles.append(validated_allele)
            assert (
                validated_allele.gene == mhc1.name.name
            ), "The allele referring to gene {} is inside gene {}".format(
                validated_allele.gene, mhc1.name.name
            )
        mhc1.alleles = validated_alleles
        return mhc1

    @staticmethod
    def _validate_mhc2(mhc2: Mhc2) -> Mhc2:
        assert mhc2.name in Mhc2Name, "Invalid MHC II name"
        genes = mhc2.genes
        for gene in genes:
            assert gene.name in Mhc2GeneName, "Invalid gene name from MHC II"
            assert gene.name in GENES_BY_MOLECULE.get(
                mhc2.name
            ), "Gene {} referring to isoform {}".format(gene.name, mhc2.name.name)
            assert gene.zygosity in Zygosity, "Invalid zygosity"
            alleles = gene.alleles
            if gene.zygosity == Zygosity.HOMOZYGOUS:
                assert (
                    len(alleles) == 1
                ), "An homozygous gene must have 1 allele and not {}".format(
                    len(alleles)
                )
            elif gene.zygosity in [Zygosity.HETEROZYGOUS, Zygosity.HEMIZYGOUS]:
                assert (
                    len(alleles) == 2
                ), "An heterozygous or hemizygous gene must have 2 alleles and not {}".format(
                    len(alleles)
                )
            elif gene.zygosity == Zygosity.LOSS:
                assert (
                    len(alleles) == 0
                ), "An lost gene must have 0 alleles and not {}".format(len(alleles))
            validated_alleles = []
            for allele in alleles:
                validated_allele = ModelValidator.validate_mhc_allele_representation(
                    allele
                )
                validated_alleles.append(validated_allele)
                assert (
                    validated_allele.gene == gene.name.name
                ), "The allele referring to gene {} is inside gene {}".format(
                    validated_allele.gene, gene.name.name
                )
            gene.alleles = validated_alleles
        isoforms = mhc2.isoforms
        validated_isoforms = []
        for isoform in isoforms:
            validated_isoform = ModelValidator.validate_mhc2_isoform_representation(
                isoform
            )
            validated_isoforms.append(validated_isoform)
            if mhc2.name != Mhc2Name.DR:
                assert validated_isoform.alpha_chain.name in [
                    a.name for g in genes for a in g.alleles
                ], "Alpha chain allele not present in th list of alleles"
            assert validated_isoform.beta_chain.name in [
                a.name for g in genes for a in g.alleles
            ], "Beta chain allele not present in th list of alleles"
        mhc2.isoforms = validated_isoforms
        return mhc2

    @staticmethod
    def validate_mhc_allele_representation(allele: MhcAllele) -> MhcAllele:
        try:
            full_name = None
            if allele.full_name:
                # infers gene, group and protein from the name
                match = HLA_ALLELE_PATTERN.match(allele.full_name)
                assert (
                    match is not None
                ), "Allele does not match HLA allele pattern {}".format(allele.name)
                gene = match.group(1)
                group = match.group(2)
                protein = match.group(3)
                full_name = allele.full_name
            elif allele.name:
                # infers gene, group and protein from the name
                match = HLA_ALLELE_PATTERN.match(allele.name)
                assert (
                    match is not None
                ), "Allele does not match HLA allele pattern {}".format(allele.name)
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
                    "HLA allele missing required fields, either name or gene, group and protein must be provided"
                )

            assert gene in list(Mhc1Name.__members__.keys()) + list(
                Mhc2GeneName.__members__.keys()
            ), "Gene not from classic MHC: {}".format(gene)
            # builds the final allele representation and validates it just in case
            name = "HLA-{gene}*{serotype}:{protein}".format(
                gene=gene, serotype=group, protein=protein
            )
            match = HLA_ALLELE_PATTERN.match(name)
            assert (
                match is not None
            ), "Allele does not match HLA allele pattern {}".format(name)
        except AssertionError as e:
            raise NeofoxDataValidationException(e)

        return MhcAllele(
            full_name=full_name if full_name else name,
            name=name,
            gene=gene,
            group=group,
            protein=protein,
        )

    @staticmethod
    def validate_mhc2_isoform_representation(isoform: Mhc2Isoform) -> Mhc2Isoform:
        try:
            if isoform.name:
                # infers alpha and beta chains
                match = HLA_MOLECULE_PATTERN.match(isoform.name)
                if match:
                    alpha_chain = ModelValidator.validate_mhc_allele_representation(
                        MhcAllele(name=match.group(1))
                    )
                    beta_chain = ModelValidator.validate_mhc_allele_representation(
                        MhcAllele(name=match.group(2))
                    )
                else:
                    match = HLA_DR_MOLECULE_PATTERN.match(isoform.name)
                    assert (
                        match is not None
                    ), "Molecule does not match HLA isoform pattern {}".format(
                        isoform.name
                    )
                    alpha_chain = MhcAllele()
                    beta_chain = ModelValidator.validate_mhc_allele_representation(
                        MhcAllele(name=match.group(1))
                    )
            elif isoform.alpha_chain and isoform.beta_chain:
                # infers name from gene, group and protein
                alpha_chain = ModelValidator.validate_mhc_allele_representation(
                    isoform.alpha_chain
                )
                beta_chain = ModelValidator.validate_mhc_allele_representation(
                    isoform.beta_chain
                )
            else:
                raise NeofoxDataValidationException(
                    "HLA isoform missing required fields"
                )

            # builds the final allele representation and validates it just in case
            name = get_mhc2_isoform_name(alpha_chain, beta_chain)
            match = HLA_MOLECULE_PATTERN.match(name)
            match2 = HLA_DR_MOLECULE_PATTERN.match(name)
            assert (
                match is not None or match2 is not None
            ), "Molecule does not match HLA isoform pattern {}".format(name)
        except AssertionError as e:
            raise NeofoxDataValidationException(e)

        return Mhc2Isoform(name=name, alpha_chain=alpha_chain, beta_chain=beta_chain)

    @staticmethod
    def _validate_expression_values(neoantigen):
        assert (
            neoantigen.rna_expression is None or neoantigen.rna_expression >= 0
        ), "RNA expression should be a positive integer or zero {}".format(
            neoantigen.rna_expression
        )
        ModelValidator._validate_vaf(neoantigen.dna_variant_allele_frequency)
        ModelValidator._validate_vaf(neoantigen.rna_variant_allele_frequency)

    @staticmethod
    def _validate_mutation(mutation: Mutation) -> Mutation:
        assert mutation.mutated_xmer is not None and len(mutation.mutated_xmer) > 0, \
            "Missing mutated xmer on mutation"
        mutation.mutated_xmer = "".join(
            [ModelValidator._validate_aminoacid(aa) for aa in mutation.mutated_xmer]
        )
        # avoids this validation when there is no wild type
        if mutation.wild_type_xmer:
            mutation.wild_type_xmer = "".join(
                [
                    ModelValidator._validate_aminoacid(aa)
                    for aa in mutation.wild_type_xmer
                ]
            )
            mutation.position = EpitopeHelper.mut_position_xmer_seq(mutation=mutation)
        return mutation

    @staticmethod
    def _validate_vaf(vaf):
        assert (
            vaf is None or vaf == -1.0 or 0.0 <= vaf <= 1.0
        ), "VAF should be a positive integer or zero {}".format(vaf)

    @staticmethod
    def _validate_aminoacid(aminoacid):
        assert aminoacid is not None, "Amino acid field cannot be empty"
        aminoacid = aminoacid.strip()
        assert isinstance(aminoacid, str), "Amino acid has to be a string"
        # this chunk is unused but let's leave in case it is handy in the future
        if len(aminoacid) == 3:
            assert (
                aminoacid in IUPACData.protein_letters_3to1_extended.keys()
            ), "Non existing 3 letter amino acid {}".format(aminoacid)
            aminoacid = IUPACData.protein_letters_3to1_extended.get(aminoacid)
        if len(aminoacid) == 1:
            aminoacid = aminoacid.upper()
            assert (
                aminoacid in ExtendedIUPACProtein.letters
            ), "Non existing aminoacid {}".format(aminoacid)
        else:
            assert False, "Invalid aminoacid {}".format(aminoacid)
        return aminoacid

    @staticmethod
    def generate_neoantigen_identifier(neoantigen: Neoantigen) -> str:
        neoantigen.identifier = None  # this needs to be done otherwise we cannot rreproduce the id after it is set
        return base64.b64encode(
            hashlib.md5(neoantigen.to_json().encode("utf8")).digest()
        ).decode("utf8")
