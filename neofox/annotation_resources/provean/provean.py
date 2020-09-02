import re
import pysam
import gzip
from logzero import logger


class ProveanAnnotator(object):

    def __init__(self, provean_file):
        """
        The provean file is tab separated values file compressed with bgzip and tabix indexed on the ucsc id and the
        position.
        The header of the file is as follows:
        #protein_id     position        A       C       D       E       F       G       H       I       K       L
        M       N       P       Q       R       S       T       V       W       Y       Del
        UCSC_ID UCSC_ID_stripped  UCSC_ID_POS
        """
        self.aminoacid_indexes = self._load_aminoacid_indexes(provean_file)
        self.provean = pysam.TabixFile(provean_file)

    def _load_aminoacid_indexes(self, provean_file):
        with gzip.open(provean_file, "rt") as p:
            header = p.readline().split("\t")
        return {x:i for i, x in enumerate(header)}

    def get_provean_annotation(self, protein_id, position, mutated_aminoacid):
        """
        Returns the PROVEAN score of particular mutation in a protein
        :param protein_id: ucsc protein id
        :type protein_id: str
        :param position: the position in the protein
        :type position: int
        :param mutated_aminoacid: the mutated aminoacid
        :type mutated_aminoacid: str
        :return: the provean score
        :rtype str
        """
        provean_score = None
        logger.info("Fetching the PROVEAN score at {}:{}:{}".format(protein_id, position, mutated_aminoacid))
        try:
            results = self.provean.fetch(protein_id, position - 1, position)
            provean_entry = next(results).split("\t")
            provean_score = provean_entry[self.aminoacid_indexes.get(mutated_aminoacid)]
        except (IndexError, TypeError, ValueError) as ex:
            logger.error("Error fetching the PROVEAN score at {}:{}:{}".format(protein_id, position, mutated_aminoacid))
        logger.info("Fetched a PROVEAN score of {}".format(provean_score))
        return provean_score
