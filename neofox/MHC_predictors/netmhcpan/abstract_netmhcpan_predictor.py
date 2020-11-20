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


@dataclass
class PredictedEpitope:
    "this is a common data class for both netmhcpan and netmhc2pan"
    pos: int
    hla: Union[str, Mhc2Isoform]    # for MHC I a str is enough, but for MCH II we need a complex object
    peptide: str
    affinity_score: float
    rank: float


class AbstractNetMhcPanPredictor:

    @staticmethod
    def select_best_by_rank(predictions: List[PredictedEpitope]) -> PredictedEpitope:
        """reports best predicted epitope (over all alleles). indicate by rank = true if rank score should be used.
        if rank = False, Aff(nM) is used
        In case of a tie, it chooses the first peptide in alphabetical order
        """
        return min(predictions, key=lambda p: (p.rank, p.peptide))

    @staticmethod
    def select_best_by_affinity(predictions: List[PredictedEpitope]) -> PredictedEpitope:
        """reports best predicted epitope (over all alleles). indicate by rank = true if rank score should be used.
        if rank = False, Aff(nM) is used
        In case of a tie, it chooses the first peptide in alphabetical order
        """
        return min(predictions, key=lambda p: (p.affinity_score, p.peptide))

    def filter_for_WT_epitope_position(
            self, predictions: List[
                PredictedEpitope], sequence_mut, position_mutation_epitope) -> PredictedEpitope:
        """returns wt epitope info for given mutated sequence. best wt that is allowed to bind to any allele of patient
        """
        epitopes_wt = list(filter(
           lambda p: len(p.peptide) == len(sequence_mut) and p.pos == position_mutation_epitope, predictions))
        return self.select_best_by_rank(epitopes_wt)

    def filter_binding_predictions(
            self, position_of_mutation, predictions: List[PredictedEpitope]) -> List[PredictedEpitope]:
        """filters prediction file for predicted epitopes that cover mutations"""
        return list(filter(
           lambda p: EpitopeHelper.epitope_covers_mutation(position_of_mutation, p.pos, len(p.peptide)), predictions))

    def filter_for_9mers(self, predictions: List[PredictedEpitope]) -> List[PredictedEpitope]:
        """returns only predicted 9mers"""
        return list(filter(lambda p: len(p.peptide) == 9, predictions))
