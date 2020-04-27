#!/usr/bin/env python

from logzero import logger


def add_ucsc_id_to_dict(props):
    subitution = props["substitution"]
    ucsc_id = props["UCSC_transcript"]
    try:
        ucsc_epi = ucsc_id.split(".")[0]
        pos_prot = str(''.join(i for i in subitution if i.isdigit()))
        return "_".join([ucsc_epi, str(pos_prot)])
    except ValueError:
        return "_".join([ucsc_epi, "Del"])


def add_ucsc_id_to_list(ucsc_id, subst):
    try:
        ucsc_epi = ucsc_id.split(".")[0]
        pos_prot = str(''.join(i for i in subst if i.isdigit()))
        return "_".join([ucsc_epi, str(pos_prot)])
    except ValueError:
        return "_".join([ucsc_epi, "Del"])


class ProveanAnnotator(object):

    def __init__(self, provean_file, epitope_ids):
        """
        Loads provean scores as dictionary, but only for ucsc ids that are in epitope list
        """
        logger.info("Starting load of PROVEAN matrix" + provean_file)
        self.provean_matrix = {}
        with open(provean_file) as f:
            next(f)  # skips header
            for line in f:
                w = line.rstrip().split(";")
                ucsc_id_pos = w[-1]
                if ucsc_id_pos in epitope_ids:
                    self.provean_matrix[ucsc_id_pos] = w
        logger.info("PROVEAN matrix loaded")

    def add_provean_score_from_matrix(self, props):
        '''
        This function maps Provean score on given position and for specific SNV onto epitope data set (which is in form of tuple --> header + dict of ucsc_pos_id: df row)
        '''
        aa_mut = props["MUT_AA"]
        ucsc_pos_epi = props["UCSC_ID_position"]
        head_prov = ["protein_id", "position", "A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R",
                     "S", "T", "V", "W", "Y", "Del"]
        try:
            return self.provean_matrix[ucsc_pos_epi][head_prov.index(aa_mut)]
        except (ValueError, KeyError) as e:
            return "NA"


"""
def add_provean_score(props, file_prov):
    '''
    This function maps Provean score on given position and for specific SNV onto epitope data set (which is in form of tuple --> header + dict of ucsc_pos_id: df row)
    VERY SLOW!
    '''
    with open(file_prov) as f:
        aa_mut = props["MUT_AA"]
        ucsc_pos_epi = props["UCSC_ID_position"]
        ucsc_epi = props["UCSC_transcript"]
        z = "NA"
        header = next(f)
        head_prov = ["protein_id", "position", "A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y", "Del"]
        for line in f:
            w = line.rstrip().split(";")
            prot_pos = w[-1]
            #print prot_pos
            if prot_pos  == ucsc_pos_epi:
                #print aa_mut
                try:
                    z = w[head_prov.index(aa_mut)]
                except ValueError:
                    z = "NA"
                print >> sys.stderr, z
                return z
            else:
                z = "NA"
        print >> sys.stderr, z
        return z
"""

"""
def add_provean_score_from_pandas(props, prov_d):
    '''
    This function maps Provean score on given position and for specific SNV onto epitope data set (which is in form of tuple --> header + dict of ucsc_pos_id: df row)
    provean matrix imported as pandas dataframe. VERY SLOW SEARCH!
    '''
    aa_mut = props["MUT_AA"]
    ucsc_pos_epi = props["UCSC_ID_position"]
    #print prov_d.head()
    #print prov_d.df.head()()[1:10]
    try:
        prov_row = prov_d.loc[prov_d["UCSC_ID_POS"] == ucsc_pos_epi]
        z = prov_row.iloc[0][aa_mut]
        return str(z)
    except (KeyError, IndexError) as e:
        return "NA"
"""
