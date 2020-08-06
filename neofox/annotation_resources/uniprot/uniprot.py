from logzero import logger
from Bio import SeqIO


class Uniprot(object):

    def __init__(self, fasta_proteome):
        logger.info("Loading Uniprot...")
        self.uniprot = self._load_proteome(fasta_proteome)
        logger.info("Loaded Uniprot. {} protein sequences".format(len(self.uniprot)))

    @staticmethod
    def _load_proteome(fasta_proteome):
        return [record.seq for record in SeqIO.parse(fasta_proteome, "fasta")]

    def is_sequence_not_in_uniprot(self, sequence):
        # TODO: make this function return booleans please
        is_not_in_uniprot = "1"
        if any([sequence in entry for entry in self.uniprot]):
            is_not_in_uniprot = "0"
        return is_not_in_uniprot