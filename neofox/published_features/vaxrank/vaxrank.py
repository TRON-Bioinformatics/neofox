# Python version with Biopython
# . /etc/profile.d/modules.sh; module load software/python/python-2.7.9

'''
This script takes as neofox table from iCAM pipeline and calculates Literature feature of neoantigens
'''

# import modules
import math
from typing import List

from neofox.model.neoantigen import Annotation
from neofox.model.wrappers import AnnotationFactory


class VaxRank():

    def __init__(self):
        self.total_binding_score = None
        self.ranking_score = None
        self.expression_score = None

    def logistic_epitope_score(
            self,
            ic50,
            midpoint=350.0,
            width=150.0,
            ic50_cutoff=5000.0):  # TODO: add these default values into CLI as arguments
        """
        Map from IC50 values to score where 1.0 = strong binder, 0.0 = weak binder
        Default midpoint and width for logistic determined by max likelihood fit
        for data from Alessandro Sette's 1994 paper:
           "The relationship between class I binding affinity
            and immunogenicity of potential cytotoxic T cell epitopes.
        adapted from: https://github.com/openvax/vaxrank/blob/master/vaxrank/epitope_prediction.py
        """
        if ic50 >= ic50_cutoff:
            return 0.0

        rescaled = (float(ic50) - midpoint) / width
        # simplification of 1.0 - logistic(x) = logistic(-x)
        logistic = 1.0 / (1.0 + math.exp(rescaled))

        # since we're scoring IC50 values, let's normalize the output
        # so IC50 near 0.0 always returns a score of 1.0
        normalizer = 1.0 / (1.0 + math.exp(-midpoint / width))

        return logistic / normalizer

    def total_binding(self, mut_scores):
        """
        adapted from: https://github.com/openvax/vaxrank/blob/master/vaxrank/epitope_prediction.py
        sums up MHC binding scores of all possible neoepitope candidates, transformed with logistic function into values between 0 and 1
        """
        mut_scores_logistic = []
        mut_scores_list = mut_scores.split("/")
        # print mut_scores_list

        # logistic transformation and sum over all epitopes deriving from mutations
        [mut_scores_logistic.append(self.logistic_epitope_score(ic50=float(mhc_affinity))) for mhc_affinity in
         mut_scores_list]
        # print mut_scores_logistic
        return sum(mut_scores_logistic)

    def combined_score(self):
        """
        adapted from: https://github.com/openvax/vaxrank/blob/master/vaxrank/epitope_prediction.py
        final ranking score implemented in VaxRank
        """
        # print "rank score: " + str(float(self.expression_score) * float(self.total_binding_score))
        combined_score = None
        try:
            combined_score = self.expression_score * self.total_binding_score
        except (ValueError, TypeError):
            pass
        return combined_score

    def run(self, mutation_scores, expression_score):
        self.expression_score = expression_score
        self.total_binding_score = self.total_binding(mutation_scores)
        self.ranking_score = self.combined_score()

    def get_annotations(self) -> List[Annotation]:
        return [
            AnnotationFactory.build_annotation(value=self.total_binding_score, name="vaxrank_binding_score"),
            AnnotationFactory.build_annotation(value=self.ranking_score, name="vaxrank_total_score")
        ]