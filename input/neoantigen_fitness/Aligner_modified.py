'''
Created on Jul 26, 2017

@author: Marta Luksza, mluksza@ias.edu
'''
from math import log, exp

from Bio import pairwise2
from Bio.Blast import NCBIXML
from Bio.SubsMat import MatrixInfo as matlist


class Aligner(object):
    '''
    Class to align neoantigens with IEDB epitopes and compute TCR-recognition
    probabilities.
    '''
    INF = float("inf")

    @staticmethod
    def align(seq1, seq2):
        '''
        Smith-Waterman alignment with default parameters.
        '''
        matrix = matlist.blosum62
        gap_open = -11
        gap_extend = -1
        aln = pairwise2.align.localds(seq1.upper(), seq2.upper(), matrix, gap_open, gap_extend)
        return aln

    @staticmethod
    def logSum(v):
        '''
        compute the logarithm of a sum of exponentials
        '''
        if len(v) == 0:
            return -Aligner.INF
        ma = max(v)
        if ma == -Aligner.INF:
            return -Aligner.INF
        return log(sum([exp(x - ma) for x in v])) + ma

    def __init__(self):
        # dictionary of computed Ri-values mapped to neoantigen identifiers
        self.Ri = {}
        # dictionary of IEDB epitope alignments mapped to neoantigen identifiers
        self.alignments = {}
        # dictionary of the highest scoring alignments mapped to neoantigen identifiers
        self.maximum_alignment = {}

    def readAllBlastAlignments(self, xmlpath):
        '''
        Read precomputed blastp alignments from xml files,
        compute alignment scores,
        find the highest scoring alignment for each neoantigen.
        '''
        f = open(xmlpath)
        blast_records = NCBIXML.parse(f)
        maxscore = {}
        try:
            for brecord in blast_records:
                nid = int(str(brecord.query).split("_")[1])
                for alignment in brecord.alignments:
                    if not nid in self.alignments:
                        self.alignments[nid] = {}
                        self.maximum_alignment[nid] = None
                        self.maximum_alignment[nid] = 0
                        maxscore[nid] = 0
                    species = " ".join((str(alignment).split())[1:-3])
                    for hsp in alignment.hsps:
                        if not "-" in hsp.query and not "-" in hsp.sbjct:
                            al = Aligner.align(hsp.query, hsp.sbjct)
                            if len(al) > 0:
                                al = al[0]
                                self.alignments[nid][species] = al
                                if al[2] > maxscore[nid]:
                                    self.maximum_alignment[nid] = species
                                    maxscore[nid] = al[2]
        except ValueError as e:
            pass
        f.close()

    def computeR(self, a=26, k=4.87):
        '''
        Compute TCR-recognition probabilities for each neoantigen.
        '''
        # iterate over all neoantigens
        for i in self.alignments:
            # energies of all bound states of neoantigen i
            bindingEnergies = [-k * (a - el[2]) for el in list(self.alignments[i].values())]
            # partition function, over all bound states and an unbound state
            lZ = Aligner.logSum(bindingEnergies + [0])
            lGb = Aligner.logSum(bindingEnergies)
            R = exp(lGb - lZ)
            self.Ri[i] = R

    def getR(self, i):
        '''
        Return precomputed R value and the highest scoring alignment
        for a given neoantigen i.
        '''
        emptyAlignment = [None, None, 0]
        if i in self.Ri:
            species = self.maximum_alignment[i]
            al = self.alignments[i][species]
            species = str(species).replace(" ", "_")
            return [self.Ri[i], species, al]
        return [0., None, emptyAlignment]


import sys

if __name__ == "__main__":
    a = Aligner()
    a.readAllBlastAlignments(sys.argv[1])
    a.computeR()
    n = []
    with open(sys.argv[2]) as f:
        for line in f:
            if line.startswith(">"):
                nid = int(line.strip("\n").split("_")[1])
                continue
            n.append((nid, line.strip("\n")))
    for nid, i in n:
        x = a.getR(nid)
        if x[1] != None:
            print(x)
