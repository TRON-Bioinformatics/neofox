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
from neofox import NOT_AVAILABLE_VALUE, MHC_II, MHC_I
from collections import defaultdict
import orjson as json
import numpy as np
from neofox.model.mhc_parser import MhcParser
from neofox.model.neoantigen import (
    Neoantigen,
    Patient,
    PredictedEpitope,
)
from neofox.model.factories import PatientFactory, NeoantigenFactory
from neofox.references.references import MhcDatabase


class ModelConverter(object):

    @staticmethod
    def parse_candidate_file(candidate_file: str) -> List[Neoantigen]:
        """
        :param candidate_file: the path to an neoantigen candidate input file
        :return neoantigens in model objects
        """
        data = pd.read_csv(
            candidate_file, sep="\t",
            # NOTE: forces the types of every column to avoid pandas setting the wrong type for corner cases
            dtype={
                "gene": str,
                "wildTypeXmer": str,
                "mutatedXmer": str,
                "patientIdentifier": str,
                "dnaVariantAlleleFrequency": float,
                "rnaExpression": float,
                "rnaVariantAlleleFrequency": float
            }
        )

        # NOTE: this is the support for the NeoFox format
        data = data.replace({np.nan: None})
        neoantigens = ModelConverter._neoantigens_csv2objects(data)

        return neoantigens

    @staticmethod
    def parse_candidate_neoepitopes_file(candidate_file: str, mhc_database: MhcDatabase) -> List[PredictedEpitope]:
        data = pd.read_csv(
            candidate_file, sep="\t",
            # NOTE: forces the types of every column to avoid pandas setting the wrong type for corner cases
            dtype={
                "gene": str,
                "mutatedPeptide": str,
                "wildTypePeptide": str,
                "dnaVariantAlleleFrequency": float,
                "rnaExpression": float,
                "rnaVariantAlleleFrequency": float,
                "patientIdentifier": str,
                "alleleMhcI": str,
                "isoformMhcII": str,
            }
        )

        # NOTE: this is the support for the NeoFox format
        data = data.replace({np.nan: None})
        neoepitopes = ModelConverter._neoepitopes_csv2objects(data, mhc_database)
        return neoepitopes


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
                tumor_type=patient_dict.get("tumorType"),
                mhc_alleles=patient_dict.get("mhcIAlleles", []),
                mhc2_alleles=patient_dict.get("mhcIIAlleles", []),
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
    def annotations2neoantigens_table(neoantigens: List[Neoantigen]) -> pd.DataFrame:
        dfs = []
        neoantigens_df = ModelConverter._neoantigens2table(neoantigens)
        neoantigens_df.replace({None: NOT_AVAILABLE_VALUE}, inplace=True)
        # we set the order of columns
        neoantigens_df = neoantigens_df.loc[:,
                         ["patientIdentifier",
                          "gene",
                          "mutatedXmer",
                          "wildTypeXmer",
                          "position",
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
    def annotations2epitopes_table(neoantigens: List[Neoantigen], mhc: str) -> pd.DataFrame:

        assert(mhc in [MHC_I, MHC_II], 'Bad MHC value')

        epitopes_dfs = []
        for n in neoantigens:
            # parses epitopes from a neoantigen into a data frame
            patient_identifier = n.patient_identifier
            epitopes = n.neoepitopes_mhc_i if mhc == MHC_I else n.neoepitopes_mhc_i_i
            epitopes_temp_df = ModelConverter._objects2dataframe(epitopes)
            epitopes_temp_df['patient_identifier'] = patient_identifier

            # adapts output table depending on MHC type
            if mhc == MHC_I:
                epitopes_temp_df.drop(list(epitopes_temp_df.filter(regex='isoformMhcII.*')), axis=1, inplace=True)
            else:
                epitopes_temp_df.drop(list(epitopes_temp_df.filter(regex='alleleMhcI.*')), axis=1, inplace=True)

            # annotations need a custom parsing, thus we remove these columns
            epitopes_temp_df.drop(list(epitopes_temp_df.filter(regex='neofoxAnnotations.*')), axis=1, inplace=True)

            # parses the annotations from each of the epitopes into a data frame
            annotations_dfs = []
            for e in epitopes:
                annotations = [a.to_dict() for a in e.neofox_annotations.annotations]
                annotations_temp_df = (pd.DataFrame(annotations).set_index("name").transpose())
                annotations_dfs.append(annotations_temp_df)
            if len(annotations_dfs) > 0:
                annotations_df = pd.concat(annotations_dfs, sort=True).reset_index()
                del annotations_df["index"]

                # puts together both data frames
                epitopes_temp_df = pd.concat([epitopes_temp_df, annotations_df], axis=1)

            epitopes_temp_df.replace({None: NOT_AVAILABLE_VALUE}, inplace=True)
            epitopes_dfs.append(epitopes_temp_df)

        # concatenates all together
        epitopes_df = pd.concat(epitopes_dfs)

        return epitopes_df

    @staticmethod
    def annotated_neoepitopes2epitopes_table(neoepitopes: List[PredictedEpitope], mhc: str) -> pd.DataFrame:

        assert (mhc in [MHC_I, MHC_II], 'Bad MHC value')

        epitopes_df = ModelConverter._objects2dataframe(neoepitopes)

        # adapts output table depending on MHC type
        if mhc == MHC_I:
            epitopes_df.drop(list(epitopes_df.filter(regex='isoformMhcII.*')), axis=1, inplace=True)
        else:
            epitopes_df.drop(list(epitopes_df.filter(regex='alleleMhcI.*')), axis=1, inplace=True)

        # formats annotation columns
        epitopes_df.drop(list(epitopes_df.filter(regex='neofoxAnnotations.*')), axis=1, inplace=True)
        # the position is used to pair neoepitopes coming out of netMHCpan in neoantigen mode, not of any use here
        epitopes_df.drop(["position"], axis=1, inplace=True)

        # parses the annotations from each of the epitopes into a data frame
        annotations_dfs = []
        for e in neoepitopes:
            annotations = [a.to_dict() for a in e.neofox_annotations.annotations]
            annotations_temp_df = (pd.DataFrame(annotations).set_index("name").transpose())
            annotations_dfs.append(annotations_temp_df)
        if len(annotations_dfs) > 0:
            annotations_df = pd.concat(annotations_dfs, sort=True).reset_index()
            del annotations_df["index"]

            # puts together both data frames
            epitopes_df = pd.concat([epitopes_df, annotations_df], axis=1)

        # replace None by NA
        epitopes_df.replace({None: NOT_AVAILABLE_VALUE}, inplace=True)

        return epitopes_df

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
    def _neoantigens_csv2objects(dataframe: pd.DataFrame) -> List[Neoantigen]:
        """transforms an patients CSV into a list of objects"""
        neoantigens = []
        for _, row in dataframe.iterrows():
            neoantigen_dict = row.to_dict()

            # build the external annotations from anything not from the model
            external_annotations = neoantigen_dict.copy()
            external_annotations.pop("wildTypeXmer", None)
            external_annotations.pop("mutatedXmer", None)
            external_annotations.pop("patientIdentifier", None)
            external_annotations.pop("gene", None)
            external_annotations.pop("rnaExpression", None)
            external_annotations.pop("rnaVariantAlleleFrequency", None)
            external_annotations.pop("dnaVariantAlleleFrequency", None)

            neoantigen = NeoantigenFactory.build_neoantigen(
                wild_type_xmer=neoantigen_dict.get("wildTypeXmer"),
                mutated_xmer=neoantigen_dict.get("mutatedXmer"),
                patient_identifier=neoantigen_dict.get("patientIdentifier"),
                gene=neoantigen_dict.get("gene"),
                rna_expression=neoantigen_dict.get("rnaExpression"),
                rna_variant_allele_frequency=neoantigen_dict.get("rnaVariantAlleleFrequency"),
                dna_variant_allele_frequency=neoantigen_dict.get("dnaVariantAlleleFrequency"),
                imputed_gene_expression=neoantigen_dict.get("imputedGeneExpression"),
                **external_annotations
            )
            neoantigens.append(neoantigen)

        return neoantigens

    @staticmethod
    def _neoepitopes_csv2objects(dataframe: pd.DataFrame, mhc_database: MhcDatabase) -> List[PredictedEpitope]:
        """transforms an patients CSV into a list of objects"""
        neoepitopes = []
        mhc_parser = MhcParser.get_mhc_parser(mhc_database)
        for _, row in dataframe.iterrows():
            neoepitope_dict = row.to_dict()

            # build the external annotations from anything not from the model
            external_annotations = neoepitope_dict.copy()
            external_annotations.pop("mutatedPeptide", None)
            external_annotations.pop("wildTypePeptide", None)
            external_annotations.pop("affinityMutated", None)
            external_annotations.pop("rankMutated", None)
            external_annotations.pop("affinityWildType", None)
            external_annotations.pop("rankWildType", None)
            external_annotations.pop("alleleMhcI", None)
            external_annotations.pop("alleleMhcII", None)
            external_annotations.pop("position", None)
            external_annotations.pop("patientIdentifier", None)
            external_annotations.pop("gene", None)
            external_annotations.pop("rnaExpression", None)
            external_annotations.pop("rnaVariantAlleleFrequency", None)
            external_annotations.pop("dnaVariantAlleleFrequency", None)

            mhci_allele = neoepitope_dict.get("alleleMhcI")
            mhcii_isoform = neoepitope_dict.get("isoformMhcII")
            patient_id = neoepitope_dict.get("patientIdentifier")
            if mhci_allele is not None and mhci_allele != '':
                neoepitope = PredictedEpitope(
                    mutated_peptide=neoepitope_dict.get("mutatedPeptide"),
                    wild_type_peptide=neoepitope_dict.get("wildTypePeptide"),
                    patient_identifier=patient_id,
                    allele_mhc_i=mhc_parser.parse_mhc_allele(mhci_allele),
                    gene=neoepitope_dict.get("gene"),
                    rna_expression=neoepitope_dict.get("rnaExpression"),
                    rna_variant_allele_frequency=neoepitope_dict.get("rnaVariantAlleleFrequency"),
                    dna_variant_allele_frequency=neoepitope_dict.get("dnaVariantAlleleFrequency"),
                    imputed_gene_expression=neoepitope_dict.get("imputedGeneExpression"),
                )
            elif mhcii_isoform is not None and mhcii_isoform != '':
                neoepitope = PredictedEpitope(
                    mutated_peptide=neoepitope_dict.get("mutatedPeptide"),
                    wild_type_peptide=neoepitope_dict.get("wildTypePeptide"),
                    patient_identifier=patient_id,
                    isoform_mhc_i_i=mhc_parser.parse_mhc2_isoform(mhcii_isoform),
                    gene=neoepitope_dict.get("gene"),
                    rna_expression=neoepitope_dict.get("rnaExpression"),
                    rna_variant_allele_frequency=neoepitope_dict.get("rnaVariantAlleleFrequency"),
                    dna_variant_allele_frequency=neoepitope_dict.get("dnaVariantAlleleFrequency"),
                    imputed_gene_expression=neoepitope_dict.get("imputedGeneExpression"),
                )
            elif patient_id is not None and patient_id != '':
                neoepitope = PredictedEpitope(
                    mutated_peptide=neoepitope_dict.get("mutatedPeptide"),
                    wild_type_peptide=neoepitope_dict.get("wildTypePeptide"),
                    patient_identifier=patient_id,
                    gene=neoepitope_dict.get("gene"),
                    rna_expression=neoepitope_dict.get("rnaExpression"),
                    rna_variant_allele_frequency=neoepitope_dict.get("rnaVariantAlleleFrequency"),
                    dna_variant_allele_frequency=neoepitope_dict.get("dnaVariantAlleleFrequency"),
                    imputed_gene_expression=neoepitope_dict.get("imputedGeneExpression"),
                )
            else:
                raise ValueError(
                    "Found an epitope without MHC-I allele, MHC-II isoform or patiend identifier: {}".format(
                        neoepitope_dict))
            neoepitopes.append(neoepitope)

        return neoepitopes

    @staticmethod
    def _neoantigens2table(neoantigens: List[Neoantigen]) -> pd.DataFrame:
        df = ModelConverter._objects2dataframe(neoantigens)
        df["position"] = df["position"].transform(lambda x: ",".join([str(y) for y in x]) if x is not None else x)
        return df
