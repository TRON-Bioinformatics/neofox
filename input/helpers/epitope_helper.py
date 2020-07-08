

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
        pos_mut_list = EpitopeHelper.mut_position_xmer_seq(xmer_mut=xmer_mut, xmer_wt=xmer_wt)
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
    def mut_position_xmer_seq(xmer_wt, xmer_mut):
        """
        returns position of mutation in xmer sequence. There can be more than one SNV within Xmer sequence.
        """
        pos_mut = []
        if len(xmer_wt) == len(xmer_mut):
            p1 = -1
            for i, aa in enumerate(xmer_mut):
                if aa != xmer_wt[i]:
                    p1 = i + 1
                    pos_mut.append(p1)
        else:
            p1 = 0
            # in case sequences do not have same length
            for a1, a2 in zip(xmer_wt, xmer_mut):
                if a1 == a2:
                    p1 += 1
                elif a1 != a2:
                    p1 += 1
                    pos_mut.append(p1)
        return pos_mut

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
