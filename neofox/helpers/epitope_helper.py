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

from Bio.Data import IUPACData

from neofox.helpers.blastp_runner import BlastpRunner
from neofox.model.neoantigen import PredictedEpitope, MhcAllele, Mhc2Isoform, Annotation, Annotations, \
    Neoantigen


class EpitopeHelper(object):
    @staticmethod
    def generate_nmers(neoantigen: Neoantigen, lengths, uniprot):
        """
        Generates peptides covering mutation of all lengths that are provided. Returns peptides as list
        No peptide is shorter than the minimun length provided
        There are no repetitions in the results
        """
        length_mut = len(neoantigen.mutated_xmer)
        list_peptides = set()
        for length in lengths:
            if length <= length_mut:
                starts = range(length_mut - length + 1)
                ends = [s + length for s in starts]
                for s, e in zip(starts, ends):
                    peptide = neoantigen.mutated_xmer[s:e]
                    if len(peptide) == length and uniprot.is_sequence_not_in_uniprot(peptide):
                        list_peptides.add(peptide)

        return list(list_peptides)

    @staticmethod
    def position_of_mutation_epitope(epitope: PredictedEpitope) -> int:
        """
        This function determines the position of the mutation within the epitope sequence.
        When multiple mutations are present it returns the last position
        """
        # TODO: is this efficient? No, a solution with zip is around 25% faster, maybe something else is even faster
        position = -1
        try:
            for i, aa in enumerate(epitope.mutated_peptide):
                if aa != epitope.wild_type_peptide[i]:
                    position = i + 1
        except Exception:
            position = None
        return position

    @staticmethod
    def number_of_mismatches(epitope_wild_type, epitope_mutation):
        """
        This function calculates the number of mismatches between the wt and the mutated epitope
        """
        # TODO: this is not efficient, it can be done with zip or using some library implementing Levenhstein distance
        p1 = 0
        for aa_mut, aa_wt in zip(epitope_mutation, epitope_wild_type):
            if aa_mut != aa_wt:
                p1 += 1
        return p1

    @staticmethod
    def position_in_anchor_position(position_mhci: int, peptide_length: int) -> bool:
        """
        This function determines if the mutation is located within an anchor position in mhc I.
        As an approximation, we assume that the second and the last position are anchor positions for all alleles.
        """
        anchor = None
        try:
            anchor = position_mhci == peptide_length or position_mhci == 2
        except Exception:
            pass
        return anchor

    @staticmethod
    def epitope_covers_mutation(
        position_mutation_list, position_epitope, length_epitope
    ):
        """
        checks if predicted epitope covers mutation
        """
        covers_mutation = False
        for position_mutation in position_mutation_list:
            if position_mutation != "-1":
                start = int(position_epitope)
                end = start + int(length_epitope) - 1
                if start <= position_mutation <= end:
                    covers_mutation = True
        return covers_mutation

    @staticmethod
    def contains_rare_amino_acid(peptide):
        found_rare_amino_acid = False
        for aa in peptide:
            if aa not in IUPACData.protein_letters:
                found_rare_amino_acid = True
                return found_rare_amino_acid
        return found_rare_amino_acid

    @staticmethod
    def pair_predictions(predictions, predictions_wt) -> List[PredictedEpitope]:
        for prediction in predictions:
            for prediction_wt in predictions_wt:
                if len(prediction_wt.mutated_peptide) == len(prediction.mutated_peptide) and \
                        prediction.position == prediction_wt.position and \
                        prediction.allele_mhc_i.name == prediction_wt.allele_mhc_i.name:
                    prediction.wild_type_peptide = prediction_wt.mutated_peptide
                    prediction.rank_wild_type = prediction_wt.rank_mutated
                    prediction.affinity_wild_type = prediction_wt.affinity_mutated
                    break
        return predictions

    @staticmethod
    def pair_mhcii_predictions(predictions, predictions_wt) -> List[PredictedEpitope]:
        for prediction in predictions:
            for prediction_wt in predictions_wt:
                if len(prediction_wt.mutated_peptide) == len(prediction.mutated_peptide) and \
                        prediction.position == prediction_wt.position and \
                        prediction.isoform_mhc_i_i.name == prediction_wt.isoform_mhc_i_i.name:
                    prediction.wild_type_peptide = prediction_wt.mutated_peptide
                    prediction.rank_wild_type = prediction_wt.rank_mutated
                    prediction.affinity_wild_type = prediction_wt.affinity_mutated
                    break
        return predictions

    @staticmethod
    def get_empty_epitope():
        return PredictedEpitope(
            mutated_peptide=None,
            position=None,
            allele_mhc_i=MhcAllele(name=None),
            isoform_mhc_i_i=Mhc2Isoform(name=None),
            affinity_mutated=None,
            rank_mutated=None,
        )

    @staticmethod
    def select_best_by_rank(predictions: List[PredictedEpitope]) -> PredictedEpitope:
        """
        Returns the peptide with the lowest (ie: meaning highest) rank and in case of tie first peptide on
        alphabetical order to ensure determinism
        """
        return max(predictions, key=lambda p: (-p.rank_mutated, p.mutated_peptide)) \
            if predictions is not None and len(predictions) > 0 else EpitopeHelper.get_empty_epitope()

    @staticmethod
    def select_best_by_affinity(predictions: List[PredictedEpitope], maximum=False) -> PredictedEpitope:
        """
        Returns the peptide with the highest affinity score and in case of tie first peptide on
        alphabetical order to ensure determinism
        By default the highest affinity score is the lowest (ie: netmhc family) if maximum=True then the highest
        score is the highest (ie: mixmhcpred and PRIME)
        """
        if maximum:
            return max(predictions, key=lambda p: (p.affinity_mutated, p.mutated_peptide)) \
                if predictions is not None and len(predictions) > 0 else EpitopeHelper.get_empty_epitope()
        else:
            return max(predictions, key=lambda p: (-p.affinity_mutated, p.mutated_peptide)) \
                if predictions is not None and len(predictions) > 0 else EpitopeHelper.get_empty_epitope()

    @staticmethod
    def remove_peptides_in_proteome(predictions: List[PredictedEpitope], uniprot
                                    ) -> List[PredictedEpitope]:
        """filters prediction file for predicted epitopes that cover mutations by searching for epitope
        in uniprot proteome database with an exact match search"""
        return list(
            filter(
                lambda p: uniprot.is_sequence_not_in_uniprot(
                    p.mutated_peptide
                ),
                predictions,
            )
        )

    @staticmethod
    def filter_for_9mers(predictions: List[PredictedEpitope]) -> List[PredictedEpitope]:
        """returns only predicted 9mers"""
        return list(filter(lambda p: len(p.mutated_peptide) == 9, predictions))

    @staticmethod
    def filter_peptides_covering_snv(
            position_of_mutation, predictions: List[PredictedEpitope]) -> List[PredictedEpitope]:
        """filters prediction file for predicted epitopes that cover mutations"""
        return list(
            filter(
                lambda p: EpitopeHelper.epitope_covers_mutation(
                    position_of_mutation, p.position, len(p.mutated_peptide)
                ),
                predictions,
            )
        )

    @staticmethod
    def set_wt_epitope_by_homology(predictions: List[PredictedEpitope], blastp_runner: BlastpRunner) -> List[PredictedEpitope]:
        """returns wt epitope for each neoepitope candidate of a neoantigen candidate from an alternative mutation
        class by a BLAST search."""

        for p in predictions:
            p.wild_type_peptide = blastp_runner.get_most_similar_wt_epitope(p.mutated_peptide)
        return predictions

    @staticmethod
    def get_epitope_id(epitope):
        if epitope.allele_mhc_i is not None:
            return "{}-{}".format(epitope.allele_mhc_i.name, epitope.mutated_peptide)
        elif epitope.isoform_mhc_i_i is not None:
            return "{}-{}".format(epitope.isoform_mhc_i_i.name, epitope.mutated_peptide)
        else:
            raise ValueError('Cannot build id on an epitope without HLA allele or isoform')

    @staticmethod
    def get_annotation_by_name(annotations: List[Annotation], name: str) -> str:
        result = None
        found = False
        for a in annotations:
            if a.name == name:
                result = a.value
                found = True
                break
        if not found:
            raise ValueError("Expected annotation '{}' not found".format(name))
        return result
