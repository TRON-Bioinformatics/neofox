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
from typing import List, Union
from dataclasses import dataclass
from neofox.model.neoantigen import Mhc2Isoform
from neofox.helpers.epitope_helper import EpitopeHelper
from neofox.helpers.runner import Runner
from neofox.helpers.blastp_runner import BlastpRunner
from neofox.references.references import DependenciesConfiguration
import os
from neofox.model.mhc_parser import MhcParser


@dataclass
class PredictedEpitope:
    """this is a common data class for both netmhcpan and netmhc2pan"""
    pos: int
    hla: Union[
        str, Mhc2Isoform
    ]  # for MHC I a str is enough, but for MCH II we need a complex object
    peptide: str
    affinity_score: float
    rank: float


class AbstractNetMhcPanPredictor(BlastpRunner):
    def __init__(self, runner: Runner, configuration: DependenciesConfiguration, proteome_db, mhc_parser: MhcParser):
        super().__init__(runner=runner, configuration=configuration)
        self.runner = runner
        self.configuration = configuration
        self.proteome_db = proteome_db
        self.mhc_parser = mhc_parser

    @staticmethod
    def select_best_by_rank(predictions: List[PredictedEpitope]) -> PredictedEpitope:
        """reports best predicted epitope (over all alleles). indicate by rank = true if rank score should be used.
        if rank = False, Aff(nM) is used
        In case of a tie, it chooses the first peptide in alphabetical order
        """
        return min(predictions, key=lambda p: (p.rank, p.peptide)) \
            if predictions is not None and len(predictions) > 0 else None

    @staticmethod
    def select_best_by_affinity(
        predictions: List[PredictedEpitope],
    ) -> PredictedEpitope:
        """reports best predicted epitope (over all alleles). indicate by rank = true if rank score should be used.
        if rank = False, Aff(nM) is used
        In case of a tie, it chooses the first peptide in alphabetical order
        """
        best_peptide = min(predictions, key=lambda p: (p.affinity_score, p.peptide))
        return best_peptide\
            if predictions is not None and len(predictions) > 0 else None

    @staticmethod
    def filter_wt_predictions_from_best_mutated(
        predictions: List[PredictedEpitope], mutated_prediction: PredictedEpitope
    ) -> List[PredictedEpitope]:
        """returns wt epitope info for given mutated sequence. best wt is restricted to the allele of best neoepitope"""
        return list(
            filter(
                lambda p: len(p.peptide) == len(mutated_prediction.peptide)
                and p.pos == mutated_prediction.pos
                and p.hla == mutated_prediction.hla,
                predictions,
            )
        )

    def find_wt_epitope_for_alternative_mutated_epitope(
            self,
            mutated_predictions: List[PredictedEpitope]
    ) -> List:
        """returns wt epitope for each neoepitope candidate of a neoantigen candidate from an alternative mutation
        class by a BLAST search."""
        mut_peptides = list(set([p.peptide for p in mutated_predictions]))
        wt_peptides = []
        database = os.path.join(self.proteome_db, "homo_sapiens")
        for p in mut_peptides:
            wt_peptides.append(self.get_most_similar_wt_epitope(p, database))
        wt_peptides_full = []
        for pep in mutated_predictions:
            for wt, mut in zip(wt_peptides, mut_peptides):
                if pep.peptide == mut:
                    wt_peptides_full.append(wt)
        return wt_peptides_full


    def filter_wt_predictions_from_best_mutated_alernative(
            self, mut_predictions: List[PredictedEpitope], wt_predictions: List[PredictedEpitope], best_mutated_epitope,
    ) -> PredictedEpitope:
        """returns wt epitope info for given mutated sequence. best wt is restricted to the allele of best neoepitope"""
        best_wt = None
        # TODO: do here re-parse
        for mut, wt in zip(mut_predictions, wt_predictions):
            if mut.peptide == best_mutated_epitope and mut.hla == wt.hla:
                best_wt = wt
        return best_wt

    @staticmethod
    def filter_binding_predictions(
        position_of_mutation, predictions: List[PredictedEpitope], uniprot
    ) -> List[PredictedEpitope]:
        """filters prediction file for predicted epitopes that cover mutations by searching for epitope
        in uniprot proteome database with an exact match search"""
        return list(
            filter(
                lambda p: uniprot.is_sequence_not_in_uniprot(
                    p.peptide
                ),
                predictions,
            )
        )

    def filter_for_9mers(
        self, predictions: List[PredictedEpitope]
    ) -> List[PredictedEpitope]:
        """returns only predicted 9mers"""
        return list(filter(lambda p: len(p.peptide) == 9, predictions))

    @staticmethod
    def filter_binding_predictions_wt_snv(
            position_of_mutation, predictions: List[PredictedEpitope]
    ) -> List[PredictedEpitope]:
        """filters prediction file for predicted epitopes that cover mutations"""
        return list(
            filter(
                lambda p: EpitopeHelper.epitope_covers_mutation(
                    position_of_mutation, p.pos, len(p.peptide)
                ),
                predictions,
            )
        )


