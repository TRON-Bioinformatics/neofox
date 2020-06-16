# Python version with Biopython
# . /etc/profile.d/modules.sh; module load software/python/python-2.7.9

'''
This script takes as input table from iCAM pipeline and calculates Literature feature of neoantigens
'''

# import modules
import math


class VaxRank():
    def __init__(self):
        self.total_binding_score = ""
        self.ranking_score = ""
        self.expression_score = ""

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
        total_score = sum(mut_scores_logistic)

        return str(total_score)

    def combined_score(self):
        """
        adapted from: https://github.com/openvax/vaxrank/blob/master/vaxrank/epitope_prediction.py
        final ranking score implemented in VaxRank
        """
        # print "rank score: " + str(float(self.expression_score) * float(self.total_binding_score))
        try:
            return str(float(self.expression_score) * float(self.total_binding_score))
        except ValueError:
            return "NA"

    def main(self, mutation_scores, expression_score):
        self.expression_score = expression_score
        self.total_binding_score = self.total_binding(mutation_scores)
        self.ranking_score = self.combined_score()


# if __name__ == '__main__':
#     # import epitope
#     # import predict_all_epitopes
#     # from helpers import data_import
#
#     epi_dict = {
#         "MB_affinities": "2568.0/6085.8/582.9/18868.9/26941.9/3681.9/29802.4/38322.6/26393.7/527.4/15690.1/189.1/15413.3/834.6/18125.5/22573.8/11247.4/36754.7/21621.7/7136.1/20981.4/2814.9/38172.8/1398.2/20769.8/35079.3/29867.9/21218.6/19437.3/35460.1/28858.4/31746.6/7889.4/40069.6/43955.6/7185.5/19266.3/41141.1/2568.0/6085.8/582.9/18868.9/26941.9/3681.9/29802.4/38322.6/26393.7/527.4/15690.1/189.1/15413.3/834.6/18125.5/22573.8/11247.4/36754.7/21621.7/7136.1/20981.4/2814.9/38172.8/1398.2/20769.8/35079.3/29867.9/21218.6/19437.3/35460.1/28858.4/31746.6/7889.4/40069.6/43955.6/7185.5/19266.3/41141.1/13582.6/33153.6/30396.7/44630.4/41746.0/38070.1/40505.1/45367.0/43947.6/39591.2/2692.9/10750.1/17879.5/39429.2/45230.3/33230.2/29433.2/44502.6/37430.5/42170.9/39078.0/2488.5/34630.5/43600.9/40167.3/28268.1/43753.0/32705.8/36653.8/44132.0/26697.3/44069.0/42608.9/43093.9/45060.8/35449.8/42832.2/45441.2/13582.6/33153.6/30396.7/44630.4/41746.0/38070.1/40505.1/45367.0/43947.6/39591.2/2692.9/10750.1/17879.5/39429.2/45230.3/33230.2/29433.2/44502.6/37430.5/42170.9/39078.0/2488.5/34630.5/43600.9/40167.3/28268.1/43753.0/32705.8/36653.8/44132.0/26697.3/44069.0/42608.9/43093.9/45060.8/35449.8/42832.2/45441.2/13477.5/34902.1/35727.0/44343.5/34096.2/20217.5/42885.5/44374.7/43128.4/35133.2/23984.1/14154.1/20684.8/34442.9/43631.5/9376.9/21796.7/45814.9/30521.6/39908.7/32105.2/4218.7/18279.2/37758.8/20711.9/11696.2/45202.4/27644.4/27896.2/39722.5/23205.3/41291.8/27966.3/42898.0/42222.0/19525.7/40617.0/43220.4/13477.5/34902.1/35727.0/44343.5/34096.2/20217.5/42885.5/44374.7/43128.4/35133.2/23984.1/14154.1/20684.8/34442.9/43631.5/9376.9/21796.7/45814.9/30521.6/39908.7/32105.2/4218.7/18279.2/37758.8/20711.9/11696.2/45202.4/27644.4/27896.2/39722.5/23205.3/41291.8/27966.3/42898.0/42222.0/19525.7/40617.0/43220.4",
#         "Expression_Mutated_Transcript": "3.20194922282"}
#
#     p = VaxRank()
#     p.main(epi_dict)
#     print("expression: " + p.expression_score)
#     print("binding score: " + p.total_binding_score)
#     print("total score: " + p.ranking_score)
