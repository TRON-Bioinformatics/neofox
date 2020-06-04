

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
        pos_mut = int(EpitopeHelper.mut_position_xmer_seq(xmer_mut=xmer_mut, xmer_wt=xmer_wt))
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
        returns position of mutation in xmer sequence
        """
        p1 = -1
        if len(xmer_wt) == len(xmer_mut):
            p1 = -1
            for i, aa in enumerate(xmer_mut):
                if aa != xmer_wt[i]:
                    p1 = i + 1
        else:
            p1 = 0
            # in case sequences do not have same length
            for a1, a2 in zip(xmer_wt, xmer_mut):
                if a1 == a2:
                    p1 += 1
        return str(p1)

    @staticmethod
    def epitope_covers_mutation(position_mutation, position_epitope, length_epitope):
        """
        checks if predicted epitope covers mutation
        """
        cover = False
        if position_mutation != "-1":
            start = int(position_epitope)
            end = start + int(length_epitope) - 1
            if int(position_mutation) >= start and int(position_mutation) <= end:
                cover = True
        return cover

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
