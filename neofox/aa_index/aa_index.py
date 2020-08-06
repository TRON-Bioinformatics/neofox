"""
H ANDN920101
D alpha-CH chemical shifts (Andersen et al., 1992)
R PMID:1575719
A Andersen, N.H., Cao, B. and Chen, C.
T Peptide/protein structure analysis using the chemical shift index method:
  upfield alpha-CH values reveal dynamic helices and aL sites
J Biochem. and Biophys. Res. Comm. 184, 1008-1014 (1992)
C BUNA790102    0.949
I    A/L     R/K     N/M     D/F     C/P     Q/S     E/T     G/W     H/Y     I/V
    4.35    4.38    4.75    4.76    4.65    4.37    4.29    3.97    4.63    3.95
    4.17    4.36    4.52    4.66    4.44    4.50    4.35    4.70    4.60    3.95
//
"""
import os
from logzero import logger


AA_INDEX1_FILENAME = 'aaindex1'
AA_INDEX2_FILENAME = 'aaindex2'


class AaIndex(object):

    def __init__(self):
        aaindex1_file = os.path.join(os.path.abspath(os.path.dirname(__file__)), AA_INDEX1_FILENAME)
        aaindex2_file = os.path.join(os.path.abspath(os.path.dirname(__file__)), AA_INDEX2_FILENAME)
        logger.info("Loading AA index 1...")
        self.aaindex1 = self._parse_aaindex1(aaindex1_file)
        logger.info("Loaded {} entries...".format(len(self.aaindex1)))
        logger.info("Loading AA index 2...")
        self.aaindex2 = self._parse_aaindex2(aaindex2_file)
        logger.info("Loaded {} entries...".format(len(self.aaindex2)))

    def get_aaindex1(self):
        return self.aaindex1

    def get_aaindex2(self):
        return self.aaindex2

    def _parse_aaindex1(self, fin):
        d = {}
        with open(fin) as f:
            id = ""
            keys1 = []
            keys2 = []
            values = []
            data = False
            lb = []
            for line in f:
                l = line.strip("\n")
                w = line.split()
                lb.append(l)
                if l.startswith("H"):
                    id = w[1]
                elif l.startswith("I"):
                    data = True
                    for i in w[1:]:
                        k1, k2 = i.split("/")
                        keys1.append(k1)
                        keys2.append(k2)
                elif data and l.startswith(" "):
                    for i in w:
                        if i == "NA":
                            values.append(float('nan'))
                        else:
                            values.append(float(i))
                elif l.startswith("//"):
                    dx = {}
                    for i, j in zip(keys1 + keys2, values):
                        dx[i] = j
                    d[id] = dx
                    id = ""
                    keys1 = []
                    keys2 = []
                    values = []
                    data = False
                    lb = []
        return d


    """
    H ALTS910101
    D The PAM-120 matrix (Altschul, 1991)
    R PMID:2051488
    A Altschul, S.F.
    T Amino acid substitution matrices from an information theoretic perspective
    J J. Mol. Biol. 219, 555-565 (1991)
    M rows = ARNDCQEGHILKMFPSTWYV, cols = ARNDCQEGHILKMFPSTWYV
          3.
         -3.      6.
          0.     -1.      4.
          0.     -3.      2.      5.
         -3.     -4.     -5.     -7.      9.
         -1.      1.      0.      1.     -7.      6.
          0.     -3.      1.      3.     -7.      2.      5.
          1.     -4.      0.      0.     -5.     -3.     -1.      5.
         -3.      1.      2.      0.     -4.      3.     -1.     -4.      7.
         -1.     -2.     -2.     -3.     -3.     -3.     -3.     -4.     -4.      6.
         -3.     -4.     -4.     -5.     -7.     -2.     -4.     -5.     -3.      1.      5.
         -2.      2.      1.     -1.     -7.      0.     -1.     -3.     -2.     -2.     -4.      5.
         -2.     -1.     -3.     -4.     -6.     -1.     -4.     -4.     -4.      1.      3.      0.      8.
         -4.     -4.     -4.     -7.     -6.     -6.     -6.     -5.     -2.      0.      0.     -6.     -1.      8.
          1.     -1.     -2.     -2.     -3.      0.     -1.     -2.     -1.     -3.     -3.     -2.     -3.     -5.      6.
          1.     -1.      1.      0.     -1.     -2.     -1.      1.     -2.     -2.     -4.     -1.     -2.     -3.      1.      3.
          1.     -2.      0.     -1.     -3.     -2.     -2.     -1.     -3.      0.     -3.     -1.     -1.     -4.     -1.      2.      4.
         -7.      1.     -5.     -8.     -8.     -6.     -8.     -8.     -5.     -7.     -5.     -5.     -7.     -1.     -7.     -2.     -6.     12.
         -4.     -6.     -2.     -5.     -1.     -5.     -4.     -6.     -1.     -2.     -3.     -6.     -4.      4.     -6.     -3.     -3.     -1.      8.
          0.     -3.     -3.     -3.     -2.     -3.     -3.     -2.     -3.      3.      1.     -4.      1.     -3.     -2.     -2.      0.     -8.     -3.      5.
    //
    """


    def _parse_aaindex2(self, fin):
        d = {}
        with open(fin) as f:
            id = ""
            value_lines = []
            rows = ""
            cols = ""
            data = False
            lb = []
            firstline = False
            asym = False
            for line in f:
                l = line.strip("\n")
                w = line.split()
                lb.append(l)
                if l.startswith("H"):
                    id = w[1]
                elif l.startswith("M"):
                    data = True
                    rows = w[3].strip(",")
                    cols = w[-1]
                elif data and l.startswith(" "):
                    if not firstline:
                        firstline = True
                        if len(w) == 1:
                            asym = True
                    values = []
                    for i in w:
                        if i == "NA" or i == "-":
                            values.append(float('nan'))
                        else:
                            values.append(float(i))
                    value_lines.append(values)
                elif l.startswith("//"):
                    drows = {}
                    for ri, r in enumerate(list(rows)):
                        drows[r] = {}
                        for ci, c in enumerate(list(cols)):
                            try:
                                drows[r][c] = value_lines[ri][ci]
                            except:
                                pass
                    if asym:
                        # make symetrical
                        for r in list(rows):
                            for c in list(cols):
                                if c not in drows[r]:
                                    drows[r][c] = drows[c][r]

                    d[id] = drows
                    id = ""
                    value_lines = []
                    data = False
                    lb = []
                    rows = ""
                    cols = ""
                    firstline = False
                    asym = False
        return d


if __name__ == "__main__":
    d_aaindex1 = parse_aaindex1("aa_index/aaindex1")
    d_aaindex2 = parse_aaindex2("aa_index/aaindex2")
    print(len(list(d_aaindex1.keys())))
    print(len(list(d_aaindex2.keys())))
    print(d_aaindex2["VOGG950101"])
    print(d_aaindex2["KOSJ950101"])
    print(d_aaindex2["VOGG950101"]["A"]["C"], d_aaindex2["VOGG950101"]["C"]["A"])
    print(d_aaindex2["KOSJ950101"]["A"]["C"], d_aaindex2["KOSJ950101"]["C"]["A"])

# read trompapep
# annotate peptides with:
# z descriptors for mutant plus z describtors for WT
# matrices for substitution

# R: PCA, mark immunogenic candidates ?
#    define conserved and induced binders ?
#    retest with subgroups ?

# self similarity to proteome
#   k-mer based counting - how many k-mers of mut peptide are in proteome ?
#   maybe correalte with EPAT db ?
