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
from neofox import NOT_AVAILABLE_VALUE
from collections import defaultdict
import orjson as json
import numpy as np
from neofox.model.neoantigen import (
    Neoantigen,
    Patient,
    Annotation,
)
from neofox.model.factories import PatientFactory, NeoantigenFactory
from neofox.references.references import MhcDatabase

FIELD_VAF_DNA = "VAF_in_tumor"
FIELD_VAF_RNA = "VAF_in_RNA"
FIELD_TRANSCRIPT_EXPRESSION = "transcript_expression"
FIELD_GENE = "gene"
FIELD_WILD_TYPE_XMER = "[WT]_+-13_AA_(SNV)_/_-15_AA_to_STOP_(INDEL)"
FIELD_MUTATED_XMER = "+-13_AA_(SNV)_/_-15_AA_to_STOP_(INDEL)"


class ModelConverter(object):

    @staticmethod
    def parse_candidate_file(candidate_file: str, patient_id: str = None) -> List[Neoantigen]:
        """
        :param candidate_file: the path to an neoantigen candidate input file
        :param patient_id: the patient identifier for all neoantigens in the input file, if not provided it is
        expected as column named `patient.id` or `patient`
        :return neoantigens in model objects
        """
        data = pd.read_csv(
            candidate_file, sep="\t",
            # NOTE: forces the types of every column to avoid pandas setting the wrong type for corner cases
            dtype={
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
            # NOTE: this is the support for the iCaM format
            data = data.replace({np.nan: None})
            neoantigens = []
            for _, candidate_entry in data.iterrows():
                neoantigen = ModelConverter._candidate_entry2model(
                    candidate_entry, patient_id=patient_id
                )
                neoantigen.external_annotations = [
                    # NOTE: we need to exclude the field gene from the external annotations as it matches a field
                    # in the model and thus it causes a conflict when both are renamed to gene_x and gene_y when
                    # joining
                    Annotation(name=name, value=str(value)) for name, value
                    in candidate_entry.iteritems() if name != FIELD_GENE
                ]
                neoantigens.append(neoantigen)
        else:
            # NOTE: this is the support for the NeoFox format
            data = data.replace({np.nan: None})
            neoantigens = ModelConverter._neoantigens_csv2objects(data)
        return neoantigens

    @staticmethod
    def parse_patients_file(patients_file: str, mhc_database: MhcDatabase) -> List[Patient]:
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

        patients = []
        for _, row in df.iterrows():
            patient_dict = row.to_dict()
            patient = PatientFactory.build_patient(
                identifier=patient_dict.get("identifier"),
                is_rna_available=patient_dict.get("isRnaAvailable", False),
                tumor_type=patient_dict.get("tumorType"),
                mhc_alleles=patient_dict["mhcIAlleles"],
                mhc2_alleles=patient_dict["mhcIIAlleles"],
                mhc_database=mhc_database
            )
            patients.append(patient)
        return patients

    @staticmethod
    def parse_neoantigens_file(neoantigens_file):
        return ModelConverter._neoantigens_csv2objects(pd.read_csv(neoantigens_file, sep="\t"))

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
    def objects2json(model_objects: List[betterproto.Message]):
        """
        :param model_objects: list of objects of subclass of betterproto.Message
        """
        return [o.to_dict(casing=Casing.SNAKE) for o in model_objects]

    @staticmethod
    def annotations2table(neoantigens: List[Neoantigen]) -> pd.DataFrame:
        dfs = []
        neoantigens_df = ModelConverter._neoantigens2table(neoantigens)
        neoantigens_df.replace({None: NOT_AVAILABLE_VALUE}, inplace=True)
        # we set the order of columns
        neoantigens_df = neoantigens_df.loc[:,
                         ["patientIdentifier",
                          "gene",
                          "mutation.mutatedXmer",
                          "mutation.wildTypeXmer",
                          "mutation.position",
                          "dnaVariantAlleleFrequency",
                          "rnaVariantAlleleFrequency",
                          "rnaExpression",
                          "imputedGeneExpression"
                          ]]
        for n in neoantigens:
            annotations = [a.to_dict() for a in n.neofox_annotations.annotations]
            annotations.extend([a.to_dict() for a in n.external_annotations])
            df = (
                pd.DataFrame(annotations)
                    .set_index("name")
                    .transpose()
            )

            dfs.append(df)
        neofox_annotations_df = pd.concat(dfs, sort=True).reset_index()
        del neofox_annotations_df["index"]
        df = pd.concat([neoantigens_df, neofox_annotations_df], axis=1)
        df.replace('None', NOT_AVAILABLE_VALUE, inplace=True)
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
    def _objects2dataframe(model_objects: List[betterproto.Message]) -> pd.DataFrame:
        """
        :param model_objects: list of objects of subclass of betterproto.Message
        """
        return pd.json_normalize(
            data=[n.to_dict(include_default_values=True) for n in model_objects]
        )

    @staticmethod
    def _candidate_entry2model(candidate_entry: dict, patient_id: str) -> Neoantigen:
        """parses an row from a candidate file into a model object"""

        vaf_rna_raw = candidate_entry.get(FIELD_TRANSCRIPT_EXPRESSION)
        return NeoantigenFactory.build_neoantigen(
            wild_type_xmer=candidate_entry.get(FIELD_WILD_TYPE_XMER),
            mutated_xmer=candidate_entry.get(FIELD_MUTATED_XMER),
            patient_identifier=patient_id if patient_id else candidate_entry.get("patient"),
            gene=candidate_entry.get(FIELD_GENE),
            rna_expression=vaf_rna_raw if vaf_rna_raw is not None and vaf_rna_raw >= 0 else None,
            rna_variant_allele_frequency=candidate_entry.get(FIELD_VAF_RNA),
            dna_variant_allele_frequency=candidate_entry.get(FIELD_VAF_DNA)
        )

    @staticmethod
    def _neoantigens_csv2objects(dataframe: pd.DataFrame) -> List[Neoantigen]:
        """transforms an patients CSV into a list of objects"""
        neoantigens = []
        for _, row in dataframe.iterrows():
            nested_dict = ModelConverter._flat_dict2nested_dict(flat_dict=row.to_dict())

            # build the external annotations from anything not from the model
            external_annotations = nested_dict.copy()
            external_annotations.pop("mutation", None)
            external_annotations.pop("patientIdentifier", None)
            external_annotations.pop("gene", None)
            external_annotations.pop("rnaExpression", None)
            external_annotations.pop("rnaVariantAlleleFrequency", None)
            external_annotations.pop("dnaVariantAlleleFrequency", None)

            neoantigen = NeoantigenFactory.build_neoantigen(
                wild_type_xmer=nested_dict.get("mutation", {}).get("wildTypeXmer"),
                mutated_xmer=nested_dict.get("mutation", {}).get("mutatedXmer"),
                patient_identifier=nested_dict.get("patientIdentifier"),
                gene=nested_dict.get("gene"),
                rna_expression=nested_dict.get("rnaExpression"),
                rna_variant_allele_frequency=nested_dict.get("rnaVariantAlleleFrequency"),
                dna_variant_allele_frequency=nested_dict.get("dnaVariantAlleleFrequency"),
                imputed_gene_expression=nested_dict.get("imputedGeneExpression"),
                **external_annotations
            )
            neoantigens.append(neoantigen)

        return neoantigens

    @staticmethod
    def _neoantigens2table(neoantigens: List[Neoantigen]) -> pd.DataFrame:
        df = ModelConverter._objects2dataframe(neoantigens)
        df["mutation.position"] = df["mutation.position"].transform(
            lambda x: ",".join([str(y) for y in x]) if x is not None else x)
        return df

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

