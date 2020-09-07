


class AbstractNetMhcPanPredictor(object):

    @staticmethod
    def add_best_epitope_info(epitope_tuple, column_name):
        '''returns desired information of prediction of best epitope from netmhcpan output;
        e.g. "%Rank": MHC I score, "HLA": HLA allele, "Icore": best epitope
        '''
        dat_head = epitope_tuple[0]
        dat = epitope_tuple[1]
        val = dat_head.index(column_name)
        result = None
        try:
            result = dat[val]
        except IndexError:
            pass
        return result
