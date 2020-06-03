#!/usr/bin/env python


class AbstractMixMHCpred:

    @staticmethod
    def generate_nmers(xmer_wt, xmer_mut, lengths):
        """
        Generates peptides covering mutation of all lengths that are provided. Returns peptides as list
        No peptide is shorter than the minimun length provided
        There are no repetitions in the results
        """
        length_mut = len(xmer_mut)
        list_peptides = []
        pos_mut = int(AbstractMixMHCpred.mut_position_xmer_seq(xmer_mut=xmer_mut, xmer_wt=xmer_wt))
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

    def read_mixmhcpred(self, outtmp):
        """
        imports output of MixMHCpred prediction
        """
        counter = 0
        header = []
        dat = []
        with open(outtmp) as f:
            for line in f:
                line = line.rstrip().lstrip()
                if line:
                    if line.startswith("#"):
                        continue
                    if line.startswith("Peptide"):
                        counter += 1
                        header = line.split()
                        continue
                    line = line.split()
                    dat.append(line)
        return header, dat

    def add_best_epitope_info(self, epitope_tuple, column_name):
        """
        returns desired information of prediction of best epitope from netmhcpan output;
        e.g. "%Rank": MHC I score, "HLA": HLA allele, "Icore": best epitope
        """
        dat_head = epitope_tuple[0]
        dat = epitope_tuple[1]
        val = dat_head.index(column_name)
        try:
            return dat[val]
        except IndexError:
            return "NA"

    def extract_WT_for_best(self, xmer_wt, xmer_mut, best_mut_seq):
        """
        extracts the corresponding WT epitope for best predicted mutated epitope
        """
        start = xmer_mut.find(best_mut_seq)
        l = len(best_mut_seq)
        wt_epi = xmer_wt[start:(start + l)]
        return (wt_epi)

    def extract_WT_info(self, epitope_tuple, column_name):
        """
        :param epitope_tuple:
        :param column_name:
        :return:
        """
        dat_head = epitope_tuple[0]
        dat = epitope_tuple[1][0]
        val = dat_head.index(column_name)
        try:
            return dat[val]
        except IndexError:
            return "NA"
