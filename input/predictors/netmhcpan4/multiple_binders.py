import numpy as np
import scipy.stats as stats

from logzero import logger

from input import MHC_I, MHC_II


class MultipleBinding:

    def get_means(self, list_numbers):
        """
        returns list of arithmetic, harmonic and geometric mean from a list of numbers
        """
        results = ["NA", "NA", "NA"]
        if list_numbers is not None and len(list_numbers) > 0:
            # TODO: ensure that floats are parsed before calling this method so this conversion is not needed
            list_floats = [float(x) for x in list_numbers]
            results = [
                np.mean(list_floats),
                stats.hmean(list_floats),
                stats.gmean(list_floats)
            ]
        return results

    def generate_epi_tuple(self, prediction_out, mhc=MHC_I):
        """
        Takes netmhcpan4 output or netmhcpanII output as input (parsed with Netmhc[II]panBestPrediction().filter_binding_predictions) and
        returns tuple of mhc binding rank scores, epitope and HLA allele for all predicted epitopes as list
        """
        pred_data = prediction_out[1]
        list_of_tuples = []
        for ii, i in enumerate(pred_data):
            if mhc == MHC_II:
                # rank, affinity, epitope sequence, allele
                list_of_tuples.append((i[9], i[8], i[2], i[1]))
            else:
                # rank, affinity, epitope sequence, allele
                list_of_tuples.append((i[13], i[12], i[2], i[1]))
        return list_of_tuples

    def extract_top10_epis(self, tuple_epis):
        """
        this function sorts the predicted epitopes based on the mhc rank score and returns the top10 with lowest mhc binding score
        """
        tuple_epis.sort(key=lambda x: float(x[0]))
        return tuple_epis[0:9]

    def check_for_homozygosity(self, patient_alleles):
        """
        returns alleles that occur more than one time in list of patient alleles and hence are homozygous alleles. Otherwise retunrs empty list
        """
        return [allele for allele in patient_alleles if patient_alleles.count(allele) > 1]

    def extract_best_epi_per_alelle(self, tuple_epis, alleles):
        """
        this function returns the predicted epitope with the lowest binding score for each patient allele, considering homozyogosity
        """
        homo_alleles = self.check_for_homozygosity(alleles)
        logger.info(homo_alleles)
        dict_allels = {}
        for allele in alleles:
            for epi in tuple_epis:
                if allele == epi[-1]:
                    if allele not in dict_allels:
                        dict_allels[allele] = [epi]
                    else:
                        dict_allels[allele].append(epi)
        best_epis_per_allele = []
        for allele in dict_allels:
            dict_allels[allele].sort(key=lambda x: float(x[0]))
            best_epis_per_allele.append(dict_allels[allele][0])
            if allele in homo_alleles:
                # append homozygous allleles two times
                homo_numbers = homo_alleles.count(allele)
                logger.debug(allele)
                logger.debug(homo_numbers)
                if homo_numbers == 1:
                    best_epis_per_allele.append(dict_allels[allele][0])
                else:
                    homo_best_epi = dict_allels[allele][0]
                    homo_best_epi_all = []
                    # allele already one time represented in list --> add n-t times
                    [homo_best_epi_all.append(tuple(homo_best_epi)) for i in range(homo_numbers - 1)]
                    best_epis_per_allele.extend(tuple(homo_best_epi_all))
        logger.info(best_epis_per_allele)
        return best_epis_per_allele

    def scores_to_list(self, tuple_epis):
        """
        Takes list of epitope tuple as input and returns a list of mhc rank scores of these tuples
        """
        return [epi[0] for epi in tuple_epis]

    def affinities_to_list(self, tuple_epis):
        """
        Takes list of epitope tuple as input and returns a list of mhc rank scores of these tuples
        """
        return [epi[1] for epi in tuple_epis]

    def determine_number_of_binders(self, list_scores, threshold=2):
        """
        Determines the number of HLA binders per mutation based on a threshold. Default is set to 2, which is threshold for weak binding using netmhcpan4.
        """
        number_binders = 0
        for score in list_scores:
            if float(score) < threshold:
                number_binders += 1
        number_binders = number_binders if not len(list_scores) == 0 else "NA"
        return number_binders
