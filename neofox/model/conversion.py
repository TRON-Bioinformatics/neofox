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
from logzero import logger
from collections import defaultdict
import orjson as json
import numpy as np

from neofox.model.mhc_parser import MhcParser
from neofox.model.validation import ModelValidator
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
from neofox.model.wrappers import get_mhc2_isoform_name, NOT_AVAILABLE_VALUE
from neofox.exceptions import NeofoxInputParametersException
from neofox.references.references import ReferenceFolder, HlaDatabase

EXTERNAL_ANNOTATIONS_NAME = "External"
FIELD_VAF_DNA = "VAF_in_tumor"
FIELD_VAF_RNA = "VAF_in_RNA"
FIELD_TRANSCRIPT_EXPRESSION = "transcript_expression"
FIELD_GENE = "gene"
FIELD_WILD_TYPE_XMER = "[WT]_+-13_AA_(SNV)_/_-15_AA_to_STOP_(INDEL)"
FIELD_MUTATED_XMER = "+-13_AA_(SNV)_/_-15_AA_to_STOP_(INDEL)"
GENES_BY_MOLECULE = {
    Mhc2Name.DR: [Mhc2GeneName.DRB1],
    Mhc2Name.DP: [Mhc2GeneName.DPA1, Mhc2GeneName.DPB1],
    Mhc2Name.DQ: [Mhc2GeneName.DQA1, Mhc2GeneName.DQB1],
}


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
        data = pd.read_csv(candidate_file, sep="\t",
                           # NOTE: forces the types of every column to avoid pandas setting the wrong type for corner cases
                           dtype={
                               "identifier": str,
                               "gene": str,
                               "mutation.wildTypeXmer": str,
                               "mutation.mutatedXmer": str,
                               "patientIdentifier": str,
                               "dnaVariantAlleleFrequency": float,
                               "rnaExpression": float,
                               "rnaVariantAlleleFrequency": float
                           }
                           )

        # check format of input file
        if FIELD_MUTATED_XMER in data.columns.values.tolist():
            data = data.replace({np.nan: None})
            identifier = []
            for row in data.to_json(orient='records', lines=True).splitlines():
                id = ModelConverter.generate_neoantigen_identifier(row)
                identifier.append(id)
            data['identifier'] = identifier
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
        else:
            data = data.replace({np.nan: None})
            neoantigens, external_annotations = ModelConverter.parse_neoantigens_dataframe(data)

        return neoantigens, external_annotations

    @staticmethod
    def parse_patients_file(patients_file: str, hla_database: HlaDatabase) -> List[Patient]:
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
        return ModelConverter.patient_metadata_csv2objects(df, hla_database)

    @staticmethod
    def parse_neoantigens_file(neoantigens_file):
        return ModelConverter.parse_neoantigens_dataframe(pd.read_csv(neoantigens_file, sep="\t"))

    @staticmethod
    def parse_neoantigens_dataframe(
        dataframe: pd.DataFrame,
    ) -> Tuple[List[Neoantigen], List[NeoantigenAnnotations]]:
        """
        :param dataframe: a pandas data frame with neoantigens data
        :return: the parsed CSV into model objects
        """
        return ModelConverter.neoantigens_csv2objects(
            dataframe
        )

    @staticmethod
    def parse_neoantigens_json_file(neoantigens_json_file: str) -> List[Neoantigen]:
        """
        :param neoantigens_json_file: the file to neoantigens data JSON file
        :return: the parsed JSON into model objects
        """
        return [
            Neoantigen().from_dict(n) for n in json.loads(open(neoantigens_json_file).read())
        ]

    @staticmethod
    def objects2dataframe(model_objects: List[betterproto.Message]) -> pd.DataFrame:
        """
        :param model_objects: list of objects of subclass of betterproto.Message
        """
        return pd.json_normalize(
            data=[n.to_dict(include_default_values=True) for n in model_objects]
        )

    @staticmethod
    def objects2json(model_objects: List[betterproto.Message]):
        """
        :param model_objects: list of objects of subclass of betterproto.Message
        """
        return [o.to_dict(casing=Casing.SNAKE) for o in model_objects]

    @staticmethod
    def object2series(model_object: betterproto.Message) -> pd.Series:
        """
        :param model_object: object of subclass of betterproto.Message
        """
        return pd.json_normalize(
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
    def patient_metadata_csv2objects(dataframe: pd.DataFrame, hla_database: HlaDatabase) -> List[Patient]:
        """transforms an patients CSV into a list of objects"""
        patients = []
        for _, row in dataframe.iterrows():
            patient_dict = row.to_dict()
            patient = Patient().from_dict(patient_dict)
            mhc_alleles = patient_dict["mhcIAlleles"]
            # NOTE: during the parsing of empty columns empty lists become a list with one empty string ...
            if len(mhc_alleles) > 1 or (len(mhc_alleles) == 1 and len(mhc_alleles[0]) > 0):
                patient.mhc1 = ModelConverter.parse_mhc1_alleles(mhc_alleles, hla_database)
            else:
                patient.mhc1 = None
            mhc2_alleles = patient_dict["mhcIIAlleles"]
            if len(mhc2_alleles) > 1 or (len(mhc2_alleles) == 1 and len(mhc2_alleles[0]) > 0):
                patient.mhc2 = ModelConverter.parse_mhc2_alleles(mhc2_alleles, hla_database)
            else:
                patient.mhc2 = None
            patients.append(patient)
        return patients

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
        neoantigen.identifier = candidate_entry.get("identifier")
        neoantigen.gene = candidate_entry.get(FIELD_GENE)
        # missing RNA expression values are represented as -1
        logger.info(neoantigen.patient_identifier)
        vaf_rna_raw = candidate_entry.get(FIELD_TRANSCRIPT_EXPRESSION)
        neoantigen.rna_expression = vaf_rna_raw if vaf_rna_raw is not None and vaf_rna_raw >= 0 else None
        neoantigen.rna_variant_allele_frequency = candidate_entry.get(FIELD_VAF_RNA)
        neoantigen.dna_variant_allele_frequency = candidate_entry.get(FIELD_VAF_DNA)

        return ModelValidator.validate_neoantigen(neoantigen)


    @staticmethod
    def neoantigens_csv2objects(
        dataframe: pd.DataFrame,
    ) -> Tuple[List[Neoantigen], List[NeoantigenAnnotations]]:
        """transforms an patients CSV into a list of objects"""
        neoantigens = []
        external_annotations = []
        identifier = []
        for row in dataframe.to_json(orient='records', lines=True).splitlines():
            id = ModelConverter.generate_neoantigen_identifier(row)
            identifier.append(id)
        dataframe['identifier'] = identifier
        for _, row in dataframe.iterrows():
            nested_dict = ModelConverter._flat_dict2nested_dict(flat_dict=row.to_dict())
            neoantigen = ModelConverter._rescueNoneValues(nested_dict)
            validated_neoantigen = ModelValidator.validate_neoantigen(neoantigen)
            neoantigens.append(validated_neoantigen)
            neoantigen_names = set(Neoantigen.__annotations__.keys())
            external_annotation_names = dict.fromkeys(nam for nam in nested_dict.keys() if stringcase.snakecase(nam) not in neoantigen_names)
            external_annotations.append(
                NeoantigenAnnotations(
                    neoantigen_identifier=validated_neoantigen.identifier,
                    annotations=[
                        Annotation(name=name, value=nested_dict.get(name))
                        for name in external_annotation_names
                    ],
                    annotator=EXTERNAL_ANNOTATIONS_NAME,
                )
            )
        return neoantigens, external_annotations

    @staticmethod
    def _rescueNoneValues(
            neoantigen_dict: dict,
    ) -> Neoantigen:
        neoantigen = Neoantigen().from_dict(neoantigen_dict)

        if ModelConverter._requires_rescue("dnaVariantAlleleFrequency", neoantigen_dict):
            neoantigen.dna_variant_allele_frequency = None

        if ModelConverter._requires_rescue("rnaExpression", neoantigen_dict):
            neoantigen.rna_expression = None

        if ModelConverter._requires_rescue("rnaVariantAlleleFrequency", neoantigen_dict):
            neoantigen.rna_variant_allele_frequency = None

        if ModelConverter._requires_rescue("imputedGeneExpression", neoantigen_dict):
            neoantigen.imputed_gene_expression = None

        return neoantigen

    @staticmethod
    def _requires_rescue(name, neoantigen_dict):
        return name not in neoantigen_dict or neoantigen_dict[name] is None

    @staticmethod
    def annotations2short_wide_table(
        neoantigen_annotations: List[NeoantigenAnnotations],
        neoantigens: List[Neoantigen],
    ) -> pd.DataFrame:
        dfs = []
        neoantigens_df = ModelConverter.neoantigens2table(neoantigens)
        neoantigens_df.replace({None: NOT_AVAILABLE_VALUE}, inplace=True)
        # we set the order of columns
        neoantigens_df = neoantigens_df.loc[:,
                         ["identifier",
                          "patientIdentifier",
                          "gene",
                          "mutation.mutatedXmer",
                          "mutation.wildTypeXmer",
                          "mutation.position",
                          "dnaVariantAlleleFrequency",
                          "rnaVariantAlleleFrequency",
                          "rnaExpression",
                          "imputedGeneExpression"
                          ]]
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
        df = pd.merge(left=neoantigens_df, right=annotations_df, on="identifier")
        del df["identifier"]
        return df

    @staticmethod
    def neoantigens2table(neoantigens: List[Neoantigen]) -> pd.DataFrame:
        df = ModelConverter.objects2dataframe(neoantigens)
        df["mutation.position"] = df["mutation.position"].transform(
            lambda x: ",".join([str(y) for y in x]) if x is not None else x)
        return df

    @staticmethod
    def patients2table(patients: List[Patient]) -> pd.DataFrame:

        patients_dict = []
        for p in patients:
            patient_dict = p.to_dict(include_default_values=True)
            patient_dict["mhcIAlleles"] = ",".join([a.name for m in p.mhc1 for a in m.alleles])
            patient_dict["mhcIIAlleles"] = ",".join([a.name for m in p.mhc2 for g in m.genes for a in g.alleles])
            del patient_dict["mhc1"]
            del patient_dict["mhc2"]
            patients_dict.append(patient_dict)
        df = pd.json_normalize(data=patients_dict)
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
    def parse_mhc1_alleles(alleles: List[str], hla_database: HlaDatabase) -> List[Mhc1]:
        isoforms = []
        try:
            parsed_alleles = list(map(MhcParser(hla_database).parse_mhc_allele, alleles))
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
    def parse_mhc2_alleles(alleles: List[str], hla_database: HlaDatabase) -> List[Mhc2]:
        mhc2s = []
        try:
            parsed_alleles = list(map(MhcParser(hla_database).parse_mhc_allele, alleles))
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
            ), "MHC I allele is not valid {} at {}".format(a.gene, a.full_name)

    @staticmethod
    def _validate_mhc2_alleles(parsed_alleles: List[MhcAllele]):
        for a in parsed_alleles:
            assert (
                a.gene in Mhc2GeneName.__members__
            ), "MHC II allele is not valid {} at {}".format(a.gene, a.full_name) if a.full_name != "" else "Gene from MHC II allele is empty"

    @staticmethod
    def generate_neoantigen_identifier(json_obj) -> str:
        return base64.b64encode(hashlib.md5(json_obj.encode("utf8")).digest()).decode("utf8")


