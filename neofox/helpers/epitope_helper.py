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


class EpitopeHelper(object):

    @staticmethod
    def generate_nmers(xmer_wt, xmer_mut, lengths):
        """
        Generates peptides covering mutation of all lengths that are provided. Returns peptides as list
        No peptide is shorter than the minimun length provided
        There are no repetitions in the results
        """
        length_mut = len(xmer_mut)
        list_peptides = []
        pos_mut_list = EpitopeHelper.mut_position_xmer_seq(sequence_mut=xmer_mut, sequence_wt=xmer_wt)
        for pos_mut in pos_mut_list:
            for length in lengths:
                if length <= length_mut:
                    start_first = pos_mut - length
                    starts = [start_first + s for s in range(length)]
                    ends = [s + length for s in starts]
                    for s, e in zip(starts, ends):
                        list_peptides.append(xmer_mut[s:e])
        return list(set([x for x in list_peptides if not x == "" and len(x) >= min(lengths)]))

    @staticmethod
    def mut_position_xmer_seq(sequence_wt, sequence_mut):
        """
        returns position of mutation in xmer sequence. There can be more than one SNV within Xmer sequence.
        """
        # TODO: this is not efficient. A solution using zip is 25% faster. There may be other alternatives
        pos_mut = []
        if len(sequence_wt) == len(sequence_mut):
            p1 = -1
            for i, aa in enumerate(sequence_mut):
                if aa != sequence_wt[i]:
                    p1 = i + 1
                    pos_mut.append(p1)
        else:
            p1 = 0
            # in case sequences do not have same length
            for a1, a2 in zip(sequence_wt, sequence_mut):
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
        except:
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
        except:
            pass
        return anchor

    @staticmethod
    def epitope_covers_mutation(position_mutation_list, position_epitope, length_epitope):
        """
        checks if predicted epitope covers mutation
        """
        cover_list = [False]
        for position_mutation in position_mutation_list:
            if position_mutation != "-1":
                start = int(position_epitope)
                end = start + int(length_epitope) - 1
                if position_mutation >= start and position_mutation <= end:
                    cover_list.append(True)
        cover_mutation = any(cover_list)
        return cover_mutation

    @staticmethod
    def hamming_check_0_or_1(seq1, seq2):
        '''returns number of mismatches between 2 sequences
        '''
        errors = 0
        for i in range(len(seq1)):
            if seq1[i] != seq2[i]:
                errors += 1
                if errors >= 2:
                    return errors
        return errors
