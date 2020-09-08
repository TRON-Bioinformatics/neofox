from functools import reduce
from typing import List

from logzero import logger
import os

from neofox.model.neoantigen import Annotation
from neofox.model.wrappers import AnnotationFactory


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
                freq_dict[w[0]] = float(w[1])
        return freq_dict


class AminoacidFrequency(AbstractNmerFrequency):

    AA_FREQUENCIES_FILE_NAME = "20181108_AA_freq_prot.csv"

    def __init__(self):
        super().__init__(os.path.join(
            os.path.abspath(os.path.dirname(__file__)), self.AA_FREQUENCIES_FILE_NAME))

    def _get_frequency(self, aminoacid: str) -> float:
        return self.frequency_df.get(aminoacid)

    def _get_product_4mer_frequencies(self, sequence: str) -> float:
        """
        This function extracts 4 aa that are directed to TCR (pos 4 to pos 7 within epitope) and calculates the
        product of aa frequencies
        """
        freq_prod = None
        try:
            epi_4mer = sequence[3:7]
            epi_freqs = [self.frequency_df[aa] for aa in epi_4mer]
            freq_prod = reduce(lambda x, y: x * y, epi_freqs)
        except (TypeError, KeyError) as e:
            pass
        return freq_prod

    def get_annotations(self, aminoacid: str, sequence: str) -> List[Annotation]:
        return [
            AnnotationFactory.build_annotation(name="Frequency_mutated_amino_acid", value=self._get_frequency(aminoacid)),
            AnnotationFactory.build_annotation(name="Product_frequency_4mer_directed_TCR",
                                               value=self._get_product_4mer_frequencies(sequence))
        ]


class FourmerFrequency(AbstractNmerFrequency):

    FOURMER_FREQUENCIES_FILE_NAME = "20181108_4mer_freq.csv"

    def __init__(self):
        super().__init__(os.path.join(
            os.path.abspath(os.path.dirname(__file__)), self.FOURMER_FREQUENCIES_FILE_NAME))

    def _get_frequency_4mer(self, sequence:str) -> float:
        return self.frequency_df.get(sequence[3:7]) if len(sequence) >= 8 else None

    def get_annotations(self, sequence: str) -> List[Annotation]:
        return [
            AnnotationFactory.build_annotation(name="Frequency_of_4mer_directed_TCR", value=self._get_frequency_4mer(sequence))
        ]
