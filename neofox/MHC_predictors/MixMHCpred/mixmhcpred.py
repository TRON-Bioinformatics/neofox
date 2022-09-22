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
from neofox.helpers.runner import Runner
from neofox.model.mhc_parser import MhcParser
from neofox.model.neoantigen import Annotation, Mhc1, MhcAllele, PredictedEpitope, Neoantigen
from neofox.model.factories import AnnotationFactory
from neofox.helpers import intermediate_files
import pandas as pd
import os
from logzero import logger

from neofox.references.references import DependenciesConfiguration

ALLELE = "BestAllele"
RANK = "%Rank_bestAllele"
PEPTIDE = "Peptide"
SCORE = "Score_bestAllele"


class MixMHCpred:

    ANNOTATION_PREFIX = 'MixMHCpred'
    ANNOTATION_PREFIX_WT = 'MixMHCpred_WT'

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
            self.configuration.mix_mhc_pred_alleles_list, sep="\t"
        )
        return list(alleles["Allele"])

    def _get_mixmhc_allele_representation(self, mhc_alleles: List[MhcAllele]):
        return list(
            filter(
                lambda y: y in self.available_alleles,
                map(
                    lambda x: "{gene}{group}{protein}".format(gene=x.gene, group=x.group, protein=x.protein),
                    mhc_alleles)
            )
        )

    def _parse_mixmhcpred_output(self, filename: str) -> List[PredictedEpitope]:

        parsed_results = []
        try:
            results = pd.read_csv(filename, sep="\t", comment="#")
        except EmptyDataError:
            logger.error("Results from MixMHCpred are empty, something went wrong")
            results = pd.DataFrame()

        for _, row in results.iterrows():
            # when MixMHCpred returns no results it provides a row with the peptide and NAs for other fields
            # pandas reads NAs as float nan. Skip these
            if isinstance(row[ALLELE], str):
                parsed_results.append(
                    PredictedEpitope(
                        allele_mhc_i=self.mhc_parser.parse_mhc_allele(row[ALLELE]),
                        mutated_peptide=row[PEPTIDE],
                        affinity_mutated=float(row[SCORE]),
                        rank_mutated=float(row[RANK]),
                    ))
        return parsed_results

    def _mixmhcprediction(self, mhc_alleles: List[str], potential_ligand_sequences) -> List[PredictedEpitope]:
        """
        Performs MixMHCpred prediction for desired hla allele and writes result to temporary file.
        """
        outtmp = intermediate_files.create_temp_file(prefix="mixmhcpred", suffix=".txt")
        tmpfasta = intermediate_files.create_temp_fasta(
            potential_ligand_sequences, prefix="tmp_sequence_"
        )
        command = [
            self.configuration.mix_mhc_pred,
            "-a",
            ",".join(mhc_alleles),
            "-i",
            tmpfasta,
            "-o",
            outtmp,
        ]
        self.runner.run_command(
            cmd=command
        )
        results = self._parse_mixmhcpred_output(filename=outtmp)
        os.remove(outtmp)
        os.remove(tmpfasta)
        return results

    def run(self, neoantigen: Neoantigen, mhc: List[Mhc1], uniprot):
        """Wrapper for MHC binding prediction, extraction of best epitope and check if mutation is directed to TCR"""

        # TODO: get rid of this
        self.results = None

        # TODO: we may want to extend this to 8 to 14 bp (coordinate this with netMHCpan)
        potential_ligand_sequences = EpitopeHelper.generate_nmers(
            neoantigen=neoantigen, lengths=[8, 9, 10, 11], uniprot=uniprot
        )
        if len(potential_ligand_sequences) > 0:
            mhc1_alleles = self._get_mixmhc_allele_representation([a for m in mhc for a in m.alleles])
            if len(mhc1_alleles) > 0:
                self.results = self._mixmhcprediction(mhc1_alleles, potential_ligand_sequences)
            else:
                logger.warning("None of the MHC I alleles are supported by MixMHCpred")

    def run_peptide(self, peptide: str, allele: MhcAllele) -> PredictedEpitope:
        """Runs MixMHCpred on a single peptide"""
        result = None
        mhc1_alleles = self._get_mixmhc_allele_representation([allele])
        if len(mhc1_alleles) > 0 and 8 <= len(peptide) <= 14:
            results = self._mixmhcprediction(mhc1_alleles, [peptide])
            if results:
                result = results[0]
        else:
            logger.warning("None of the MHC I alleles are supported by MixMHCpred")
        return result

    def get_annotations(self) -> List[Annotation]:

        best_result = EpitopeHelper.select_best_by_affinity(predictions=self.results, maximum=True)
        return [
            AnnotationFactory.build_annotation(
                value=best_result.mutated_peptide, name="MixMHCpred_bestScore_peptide"
            ),
            AnnotationFactory.build_annotation(
                value=best_result.affinity_mutated, name="MixMHCpred_bestScore_score"
            ),
            AnnotationFactory.build_annotation(
                value=best_result.rank_mutated, name="MixMHCpred_bestScore_rank"
            ),
            AnnotationFactory.build_annotation(
                value=best_result.allele_mhc_i.name, name="MixMHCpred_bestScore_allele"
            ),
        ]
