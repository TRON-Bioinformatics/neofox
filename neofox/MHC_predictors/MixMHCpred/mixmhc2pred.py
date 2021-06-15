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
from logzero import logger
from neofox.exceptions import NeofoxCommandException
from pandas.errors import EmptyDataError

from neofox.helpers.epitope_helper import EpitopeHelper
from neofox.model.conversion import ModelConverter
from neofox.model.mhc_parser import MhcParser

from neofox.references.references import DependenciesConfiguration

from neofox.helpers.runner import Runner

from neofox.model.neoantigen import Annotation, Mhc2, Mhc2GeneName, MhcAllele, Mutation
from neofox.model.wrappers import AnnotationFactory, get_alleles_by_gene
from neofox.helpers import intermediate_files
import pandas as pd
import os
from logzero import logger

ALLELE = "BestAllele"
PEPTIDE = "Peptide"
RANK = "%Rank_best"


class MixMhc2Pred:
    def __init__(self, runner: Runner, configuration: DependenciesConfiguration, mhc_parser: MhcParser):
        self.runner = runner
        self.configuration = configuration
        self.available_alleles = self._load_available_alleles()
        self.mhc_parser = mhc_parser

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
    def _combine_dq_dp_alleles(list_alleles: List[str]):
        """returns patient HLA-DQ/HLA-DP allele combination that are relevant for MixMHC2pred"""
        # TODO: we need to clarify the formation of pairs here AA, BB, AB
        # TODO: what are these triplets?
        alleles_pairs = [
            "__".join([allele_1, allele_2])
            for allele_1 in list_alleles
            for allele_2 in list_alleles
            if allele_1 != allele_2
        ]
        alleles_triplets = [
            "__".join([allele_1, allele_2, allele_3])
            for allele_1 in list_alleles
            for allele_2 in list_alleles
            for allele_3 in list_alleles
            if allele_1 != allele_2 and allele_1 != allele_3 and allele_2 != allele_3
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
            self._get_mixmhc2_allele_representation(dpa1_alleles + dpb1_alleles)
        )
        dq_allele_combinations = self._combine_dq_dp_alleles(
            self._get_mixmhc2_allele_representation(dqa1_alleles + dqb1_alleles)
        )

        return [
            a
            for a in self._get_mixmhc2_allele_representation(drb1_alleles)
            + dq_allele_combinations
            + dp_allele_combinations
            if a in self.available_alleles
        ]

    def _mixmhc2prediction(
        self, mhc2: List[str], potential_ligand_sequences
    ) -> pd.DataFrame:
        """
        Performs MixMHC2pred prediction for desired hla allele and writes result to temporary file.
        """
        tmpfasta = intermediate_files.create_temp_fasta(
            potential_ligand_sequences, prefix="tmp_sequence_"
        )
        outtmp = intermediate_files.create_temp_file(
            prefix="mixmhc2pred", suffix=".txt"
        )
        cmd = [
            self.configuration.mix_mhc2_pred,
            "-a",
            " ".join(mhc2),
            "-i",
            tmpfasta,
            "-o",
            outtmp,
        ]
        self.runner.run_command(cmd)
        try:
            results = pd.read_csv(outtmp, sep="\t", comment="#")
        except EmptyDataError:
            message = "Results from MixMHC2pred are empty, something went wrong"
            logger.error(message)
            raise NeofoxCommandException(message)
        os.remove(outtmp)
        return results

    def run(self, mhc: List[Mhc2], mutation: Mutation, uniprot):
        """
        Runs MixMHC2pred:
        prediction for peptides of length 13 to 18 based on Suppl Fig. 6 a in Racle, J., et al., Nat. Biotech. (2019).
        Robust prediction of HLA class II epitopes by deep motif deconvolution of immunopeptidomes.
        """
        best_peptide = None
        best_rank = None
        best_allele = None
        potential_ligand_sequences = EpitopeHelper.generate_nmers(
            mutation=mutation, lengths=[13, 14, 15, 16, 17, 18], uniprot=uniprot
        )
        # filter mps shorter < 13aa
        filtered_sequences = list(
            filter(lambda x: len(x) >= 13, potential_ligand_sequences)
        )
        if len(filtered_sequences) > 0:
            mhc2_alleles = self.transform_hla_ii_alleles_for_prediction(mhc)
            if len(mhc2_alleles) > 0:
                results = self._mixmhc2prediction(mhc2_alleles, filtered_sequences)
                # get best result by minimum rank
                best_result = results[results[RANK] == results[RANK].min()]
                try:
                    best_peptide = best_result[PEPTIDE].iat[0]
                    best_rank = best_result[RANK].iat[0]
                    best_allele = self.mhc_parser.parse_mhc2_isoform(best_result[ALLELE].iat[0]).name
                except IndexError:
                    logger.info("MixMHC2pred returned no best result")
            else:
                logger.warning("None of the MHC II alleles are supported by MixMHC2pred")
        return best_peptide, best_rank, best_allele

    def get_annotations(self, mhc: List[Mhc2], mutation: Mutation, uniprot) -> List[Annotation]:
        best_peptide, best_rank, best_allele = self.run(mhc=mhc, mutation=mutation, uniprot=uniprot)
        return [
            AnnotationFactory.build_annotation(
                value=best_peptide, name="MixMHC2pred_best_peptide"
            ),
            AnnotationFactory.build_annotation(
                value=best_rank, name="MixMHC2pred_best_rank"
            ),
            AnnotationFactory.build_annotation(
                value=best_allele, name="MixMHC2pred_best_allele"
            ),
        ]
