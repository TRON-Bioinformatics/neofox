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
from neofox.helpers.epitope_helper import EpitopeHelper
from neofox.helpers.runner import Runner
from neofox.helpers.blastp_runner import BlastpRunner
from neofox.model.neoantigen import PredictedEpitope
from neofox.references.references import DependenciesConfiguration
from neofox.model.mhc_parser import MhcParser


class AbstractNetMhcPanPredictor:
    def __init__(self, runner: Runner, configuration: DependenciesConfiguration,
                 blastp_runner: BlastpRunner, mhc_parser: MhcParser):
        self.runner = runner
        self.configuration = configuration
        self.mhc_parser = mhc_parser
        self.blastp_runner = blastp_runner

    @staticmethod
    def select_best_by_rank(predictions: List[PredictedEpitope], none_value=None) -> PredictedEpitope:
        """reports best predicted epitope (over all alleles). indicate by rank = true if rank score should be used.
        if rank = False, Aff(nM) is used
        In case of a tie, it chooses the first peptide in alphabetical order
        """
        return min(predictions, key=lambda p: (p.rank, p.peptide)) \
            if predictions is not None and len(predictions) > 0 else none_value

    @staticmethod
    def select_best_by_affinity(predictions: List[PredictedEpitope], none_value=None) -> PredictedEpitope:
        """reports best predicted epitope (over all alleles). indicate by rank = true if rank score should be used.
        if rank = False, Aff(nM) is used
        In case of a tie, it chooses the first peptide in alphabetical order
        """
        return min(predictions, key=lambda p: (p.affinity_score, p.peptide)) \
            if predictions is not None and len(predictions) > 0 else none_value

    @staticmethod
    def filter_wt_predictions_from_best_mutated(
        predictions: List[PredictedEpitope], mutated_prediction: PredictedEpitope
    ) -> List[PredictedEpitope]:
        """returns wt epitope info for given mutated sequence. best wt is restricted to the allele of best neoepitope"""
        return list(
            filter(
                lambda p: mutated_prediction.peptide is not None and
                          len(p.peptide) == len(mutated_prediction.peptide) and
                          p.position == mutated_prediction.position and
                          p.hla.name == mutated_prediction.hla.name,
                predictions,
            )
        )

    def set_wt_epitope_by_homology(self, predictions: List[PredictedEpitope]) \
            -> List[PredictedEpitope]:
        """returns wt epitope for each neoepitope candidate of a neoantigen candidate from an alternative mutation
        class by a BLAST search."""
        for p in predictions:
            p.wild_type_peptide = self.blastp_runner.get_most_similar_wt_epitope(p.peptide)
        return predictions

    def filter_wt_predictions_from_best_mutated_alernative(
            self, predictions: List[PredictedEpitope], best_mutated_epitope: PredictedEpitope
    ) -> PredictedEpitope:
        """returns wt epitope info for given mutated sequence. best wt is restricted to the allele of best neoepitope"""
        best_wt = None
        for p in predictions:
            if mut.peptide == best_mutated_epitope.peptide:
                best_wt = wt
                break
        return best_wt

    @staticmethod
    def remove_peptides_in_proteome(predictions: List[PredictedEpitope], uniprot
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
    def filter_peptides_covering_snv(
            position_of_mutation, predictions: List[PredictedEpitope]
    ) -> List[PredictedEpitope]:
        """filters prediction file for predicted epitopes that cover mutations"""
        return list(
            filter(
                lambda p: EpitopeHelper.epitope_covers_mutation(
                    position_of_mutation, p.position, len(p.peptide)
                ),
                predictions,
            )
        )


