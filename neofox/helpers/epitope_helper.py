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

from neofox.model.neoantigen import Mutation


class EpitopeHelper(object):
    @staticmethod
    def generate_nmers(mutation: Mutation, lengths, uniprot):
        """
        Generates peptides covering mutation of all lengths that are provided. Returns peptides as list
        No peptide is shorter than the minimun length provided
        There are no repetitions in the results
        """
        length_mut = len(mutation.mutated_xmer)
        list_peptides = []
        for length in lengths:
            if length <= length_mut:
                starts = range(length_mut - length + 1)
                ends = [s + length for s in starts]
                for s, e in zip(starts, ends):
                    peptide = mutation.mutated_xmer[s:e]
                    if len(peptide) == length and uniprot.is_sequence_not_in_uniprot(peptide):
                        list_peptides.append(peptide)

        return list_peptides

    @staticmethod
    def mut_position_xmer_seq(mutation: Mutation) -> List[int]:
        """
        returns position (1-based) of mutation in xmer sequence. There can be more than one SNV within Xmer sequence.
        """
        # TODO: this is not efficient. A solution using zip is 25% faster. There may be other alternatives
        pos_mut = []
        if len(mutation.wild_type_xmer) == len(mutation.mutated_xmer):
            p1 = -1
            for i, aa in enumerate(mutation.mutated_xmer):
                if aa != mutation.wild_type_xmer[i]:
                    p1 = i + 1
                    pos_mut.append(p1)
        else:
            p1 = 0
            # in case sequences do not have same length
            for a1, a2 in zip(mutation.wild_type_xmer, mutation.mutated_xmer):
                if a1 == a2:
                    p1 += 1
                elif a1 != a2:
                    p1 += 1
                    pos_mut.append(p1)
        return pos_mut

    @staticmethod
    def position_of_mutation_epitope(wild_type, mutation) -> int:
        """
        This function determines the position of the mutation within the epitope sequence.
        """
        # TODO: is this efficient? No, a solution with zip is around 25% faster, maybe something else is even faster
        position = -1
        try:
            for i, aa in enumerate(mutation):
                if aa != wild_type[i]:
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
