from neofox.helpers.epitope_helper import EpitopeHelper


class AbstractMixMHCpred(EpitopeHelper):

    @staticmethod
    def read_mixmhcpred(outtmp):
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

    @staticmethod
    def add_best_epitope_info(epitope_tuple, column_name):
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
            return None

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
            return None

    def difference_score(self, mut_score, wt_score):
        """calcualated difference in MixMHCpred scores between mutated and wt"""
        result = None
        try:
            result = wt_score - mut_score
        except ValueError:
            pass
        return result
