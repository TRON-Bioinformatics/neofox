from functools import reduce
from logzero import logger
import os


class AbstractNmerFrequency(object):

    def __init__(self, frequency_file):
        logger.info("Loading aminoacid frequencies from {}...".format(frequency_file))
        self.frequency_df = self._load_nmer_frequency(frequency_file)
        logger.info("Aminoacid frequencies loaded.")

    @staticmethod
    def _load_nmer_frequency(frequency_file):
        """
        Loads file with information of frequeny of nmers
        """
        # TODO: migrate this to pandas for loading and querying
        freq_dict = {}
        with open(frequency_file) as f:
            next(f)
            for line in f:
                w = line.rstrip().split(";")
                freq_dict[w[0]] = w[1]
        return freq_dict


class AminoacidFrequency(AbstractNmerFrequency):

    AA_FREQUENCIES_FILE_NAME = "20181108_AA_freq_prot.csv"

    def __init__(self):
        super().__init__(os.path.join(
            os.path.abspath(os.path.dirname(__file__)), self.AA_FREQUENCIES_FILE_NAME))

    def get_frequency(self, aminoacid):
        return str(self.frequency_df.get(aminoacid, "NA"))

    def get_product_4mer_frequencies(self, sequence):
        '''
        This function extracts 4 aa that are directed to TCR (pos 4 to pos 7 within epitope) and calculates the product of aa frequencies
        '''
        try:
            epi_4mer = sequence[3:7]
            epi_freqs = [float(self.frequency_df[aa]) for aa in epi_4mer]
            freq_prod = reduce(lambda x, y: x * y, epi_freqs)
            return str(freq_prod)
        except (TypeError, KeyError) as e:
            return "NA"


class FourmerFrequency(AbstractNmerFrequency):

    FOURMER_FREQUENCIES_FILE_NAME = "20181108_4mer_freq.csv"

    def __init__(self):
        super().__init__(os.path.join(
            os.path.abspath(os.path.dirname(__file__)), self.FOURMER_FREQUENCIES_FILE_NAME))

    def get_frequency_4mer(self, sequence):
        return str(self.frequency_df.get(sequence[3:7], "NA")) if len(sequence) >= 8 else "NA"
