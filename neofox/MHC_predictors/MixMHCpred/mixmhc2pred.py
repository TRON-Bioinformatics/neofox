#!/usr/bin/env python
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
from pandas.errors import EmptyDataError

from neofox.helpers.epitope_helper import EpitopeHelper
from neofox.model.mhc_parser import MhcParser, get_alleles_by_gene

from neofox.references.references import DependenciesConfiguration

from neofox.helpers.runner import Runner

from neofox.model.neoantigen import Annotation, Mhc2, Mhc2GeneName, MhcAllele, PredictedEpitope, Mhc2Isoform, \
    Neoantigen
from neofox.model.factories import AnnotationFactory
from neofox.helpers import intermediate_files
import pandas as pd
import os
from logzero import logger

ALLELE = "BestAllele"
PEPTIDE = "Peptide"
RANK = "%Rank_best"


class MixMHC2pred:

    ANNOTATION_PREFIX = 'MixMHC2pred'
    ANNOTATION_PREFIX_WT = 'MixMHC2pred_WT'

    def __init__(self, runner: Runner, configuration: DependenciesConfiguration, mhc_parser: MhcParser):
        self.runner = runner
        self.configuration = configuration
        self.available_alleles = self._load_available_alleles()
        self.mhc_parser = mhc_parser

        self.results = None

    def _load_available_alleles(self):
        """
        loads file with available HLA II alllels for MixMHC2pred prediction, returns set
        :return:
        """
        alleles = pd.read_csv(
            self.configuration.mix_mhc2_pred_alleles_list, skiprows=1, sep="\t"
        )
        return list(alleles["AlleleName"])

    @staticmethod
    def _combine_dq_dp_alleles(alpha_alleles: List[str], beta_alleles: List[str]):
        """returns patient HLA-DQ/HLA-DP allele combination that are relevant for MixMHC2pred"""
        # NOTE: there are some pairs of alleles which positive/negative binding could not be deconvoluted
        # hence the triplets. In MixMHC2pred the triplets are only of the form of two alpha chains and one beta chain.
        # NOTE2: this may be gone after upgrading to MixMHC2pred
        alleles_pairs = [
            "__".join([allele_1, allele_2])
            for allele_1 in alpha_alleles
            for allele_2 in beta_alleles
        ]
        alleles_triplets = [
            "__".join([allele_1, allele_2, allele_3])
            for allele_1 in alpha_alleles
            for allele_2 in alpha_alleles
            for allele_3 in beta_alleles
            if allele_1 != allele_2
        ]
        return alleles_pairs + alleles_triplets

    @staticmethod
    def _get_mixmhc2_allele_representation(hla_alleles: List[MhcAllele]):
        return list(
            map(
                lambda x: "{gene}_{group}_{protein}".format(
                    gene=x.gene, group=x.group, protein=x.protein
                ),
                hla_alleles,
            )
        )

    @staticmethod
    def _get_mixmhc2_isoform_representation(isoform: Mhc2Isoform):

        beta_chain = MixMHC2pred._get_mixmhc2_allele_representation([isoform.beta_chain])[0]
        if isoform.alpha_chain is not None and isoform.alpha_chain.name:
            # for DR only beta chain is provided
            alpha_chain = MixMHC2pred._get_mixmhc2_allele_representation([isoform.alpha_chain])[0]
            return "{alpha}__{beta}".format(alpha=alpha_chain, beta=beta_chain)
        return beta_chain

    def transform_hla_ii_alleles_for_prediction(self, mhc: List[Mhc2]) -> List[str]:
        """
        prepares list of HLA II alleles for prediction in required format
        """
        drb1_alleles = get_alleles_by_gene(mhc, Mhc2GeneName.DRB1)
        dpa1_alleles = get_alleles_by_gene(mhc, Mhc2GeneName.DPA1)
        dpb1_alleles = get_alleles_by_gene(mhc, Mhc2GeneName.DPB1)
        dqa1_alleles = get_alleles_by_gene(mhc, Mhc2GeneName.DQA1)
        dqb1_alleles = get_alleles_by_gene(mhc, Mhc2GeneName.DQB1)

        dp_allele_combinations = self._combine_dq_dp_alleles(
            alpha_alleles=self._get_mixmhc2_allele_representation(dpa1_alleles),
            beta_alleles=self._get_mixmhc2_allele_representation(dpb1_alleles)
        )
        dq_allele_combinations = self._combine_dq_dp_alleles(
            alpha_alleles=self._get_mixmhc2_allele_representation(dqa1_alleles),
            beta_alleles=self._get_mixmhc2_allele_representation(dqb1_alleles)
        )

        return [
            a
            for a in self._get_mixmhc2_allele_representation(drb1_alleles)
            + dq_allele_combinations
            + dp_allele_combinations
            if a in self.available_alleles
        ]

    def _parse_mixmhc2pred_output(self, filename: str) -> List[PredictedEpitope]:

        parsed_results = []
        try:
            results = pd.read_csv(filename, sep="\t", comment="#")
        except EmptyDataError:
            logger.error("Results from MixMHC2pred are empty, something went wrong")
            results = pd.DataFrame()

        for _, row in results.iterrows():
            # when MixMHC2pred returns no results it provides a row with the peptide and NAs for other fields
            # pandas reads NAs as float nan. Skip these
            if isinstance(row[ALLELE], str):
                parsed_results.append(
                    PredictedEpitope(
                        isoform_mhc_i_i=self.mhc_parser.parse_mhc2_isoform(row[ALLELE]),
                        mutated_peptide=row[PEPTIDE],
                        rank_mutated=float(row[RANK]),
                        affinity_mutated=None
                    ))
        return parsed_results

    def _mixmhc2prediction(self, isoforms: List[str], potential_ligand_sequences: List[str]) -> List[PredictedEpitope]:

        tmpfasta = intermediate_files.create_temp_fasta(potential_ligand_sequences, prefix="tmp_sequence_")
        outtmp = intermediate_files.create_temp_file(prefix="mixmhc2pred", suffix=".txt")
        cmd = [
            self.configuration.mix_mhc2_pred,
            "-a",
            " ".join(isoforms),
            "-i",
            tmpfasta,
            "-o",
            outtmp,
        ]
        self.runner.run_command(cmd)
        results = self._parse_mixmhc2pred_output(filename=outtmp)
        os.remove(outtmp)
        os.remove(tmpfasta)
        return results

    def run(self, mhc: List[Mhc2], neoantigen: Neoantigen, uniprot):
        """
        Runs MixMHC2pred:
        prediction for peptides of length 13 to 18 based on Suppl Fig. 6 a in Racle, J., et al., Nat. Biotech. (2019).
        Robust prediction of HLA class II epitopes by deep motif deconvolution of immunopeptidomes.
        """
        # TODO: get rid of this
        self.results = None

        potential_ligand_sequences = EpitopeHelper.generate_nmers(
            neoantigen=neoantigen, lengths=[13, 14, 15, 16, 17, 18], uniprot=uniprot)
        # filter mps shorter < 13aa
        filtered_sequences = list(
            filter(lambda x: len(x) >= 13, potential_ligand_sequences)
        )
        if len(filtered_sequences) > 0:
            mhc2_alleles = self.transform_hla_ii_alleles_for_prediction(mhc)
            if len(mhc2_alleles) > 0:
                self.results = self._mixmhc2prediction(
                    isoforms=mhc2_alleles, potential_ligand_sequences=filtered_sequences)
            else:
                logger.warning("None of the MHC II alleles are supported by MixMHC2pred")

    def run_peptide(self, peptide: str, isoform: Mhc2Isoform) -> PredictedEpitope:
        """
        Performs MixMHC2pred prediction for desired hla allele and writes result to temporary file.
        """
        result = None
        isoform_representation = self._get_mixmhc2_isoform_representation(isoform)
        if isoform_representation in self.available_alleles:
            results = self._mixmhc2prediction(
                isoforms=[isoform_representation],
                potential_ligand_sequences=[peptide])
            if results:
                result = results[0]
        return result

    def get_annotations(self) -> List[Annotation]:
        best_result = EpitopeHelper.select_best_by_rank(predictions=self.results)
        return [
            AnnotationFactory.build_annotation(
                value=best_result.mutated_peptide, name="MixMHC2pred_bestRank_peptide"
            ),
            AnnotationFactory.build_annotation(
                value=best_result.rank_mutated, name="MixMHC2pred_bestRank_rank"
            ),
            AnnotationFactory.build_annotation(
                value=best_result.isoform_mhc_i_i.name, name="MixMHC2pred_bestRank_allele"
            ),
        ]
