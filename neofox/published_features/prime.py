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

from neofox.exceptions import NeofoxCommandException
from neofox.helpers.epitope_helper import EpitopeHelper
from neofox.helpers.runner import Runner
from neofox.model.mhc_parser import MhcParser

from neofox.model.neoantigen import Annotation, Mhc1, MhcAllele, Mutation
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


class Prime:
    def __init__(self, runner: Runner, configuration: DependenciesConfiguration, mhc_parser: MhcParser):
        self.runner = runner
        self.configuration = configuration
        self.available_alleles = self._load_available_alleles()
        self.mhc_parser = mhc_parser

    def _load_available_alleles(self):
        """
        loads file with available HLA II alllels for Prime prediction, returns set
        :return:
        """
        alleles = pd.read_csv(
            self.configuration.prime_alleles_list, sep="\t", header=None, names=['Allele']
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

    def _prime(
        self, mhc_alleles: List[str], potential_ligand_sequences
    ) -> pd.DataFrame:
        """
        Runs PRIME for desired hla allele and writes result to temporary file.
        """
        outtmp = intermediate_files.create_temp_file(prefix="prime", suffix=".txt")
        tmpfasta = intermediate_files.create_temp_fasta(
            potential_ligand_sequences, prefix="tmp_sequence_"
        )
        command = [
            self.configuration.prime,
            "-a",
            ",".join(mhc_alleles),
            "-i",
            tmpfasta,
            "-o",
            outtmp,
            "-mix",
            self.configuration.mix_mhc_pred
        ]
        self.runner.run_command(
            cmd=command
        )
        try:
            results = pd.read_csv(outtmp, sep="\t", comment="#")
        except EmptyDataError:
            message = "Results from PRIME are empty, something went wrong [{}]. MHC I alleles {}, ligands {}".format(
                " ".join(command), ",".join(mhc_alleles), potential_ligand_sequences
            )
            logger.error(message)
            results = pd.DataFrame()
        os.remove(outtmp)
        return results

    def run(self, mutation: Mutation, mhc: List[Mhc1], uniprot):
        """Wrapper PRIME prediction, extraction of best epitope per mutations"""
        best_peptide = None
        best_rank = None
        best_allele = None
        best_score = None

        if not EpitopeHelper.contains_rare_amino_acid(mutation.mutated_xmer):
            potential_ligand_sequences = EpitopeHelper.generate_nmers(
                mutation=mutation, lengths=[8, 9, 10, 11, 12, 13, 14], uniprot=uniprot
            )
            if len(potential_ligand_sequences) > 0:
                mhc1_alleles = self._get_mixmhc_allele_representation([a for m in mhc for a in m.alleles])
                if len(mhc1_alleles) > 0:
                    results = self._prime(mhc1_alleles, potential_ligand_sequences)
                    try:
                        # get best result by maximum score
                        best_result = results[results[SCORE] == results[SCORE].max()]
                        best_peptide = best_result[PEPTIDE].iat[0]
                        best_rank = best_result[RANK].iat[0]
                        # normalize the HLA allele name
                        best_allele = self.mhc_parser.parse_mhc_allele(best_result[ALLELE].iat[0]).name
                        best_score = best_result[SCORE].iat[0]
                    except (IndexError, KeyError):
                        logger.info("PRIME returned no best result")
                else:
                    logger.warning("None of the MHC I alleles are supported by PRIME")

        return best_peptide, best_rank, best_allele, best_score

    def get_annotations(self, mutation: Mutation, mhc: List[Mhc1], uniprot) -> List[Annotation]:
        best_peptide, best_rank, best_allele, best_score = self.run(
            mhc=mhc, mutation=mutation, uniprot=uniprot
        )
        return [
            AnnotationFactory.build_annotation(
                value=best_peptide, name="PRIME_best_peptide"
            ),
            AnnotationFactory.build_annotation(
                value=best_score, name="PRIME_best_score"
            ),
            AnnotationFactory.build_annotation(
                value=best_rank, name="PRIME_best_rank"
            ),
            AnnotationFactory.build_annotation(
                value=best_allele, name="PRIME_best_allele"
            ),
        ]
