from typing import List

from logzero import logger
from Bio import SeqIO

from neofox.model.neoantigen import Annotation
from neofox.model.wrappers import AnnotationFactory


class Uniprot(object):

    def __init__(self, fasta_proteome):
        logger.info("Loading Uniprot...")
        self.uniprot = self._load_proteome(fasta_proteome)
        logger.info("Loaded Uniprot. {} protein sequences".format(len(self.uniprot)))

    @staticmethod
    def _load_proteome(fasta_proteome):
        return [record.seq for record in SeqIO.parse(fasta_proteome, "fasta")]

    def is_sequence_not_in_uniprot(self, sequence) -> bool:
        is_not_in_uniprot = True
        if any([sequence in entry for entry in self.uniprot]):
            is_not_in_uniprot = False
        return is_not_in_uniprot

    def get_annotations(self, sequence_not_in_uniprot: bool) -> List[Annotation]:
        return [
            AnnotationFactory.build_annotation(name="mutation_not_found_in_proteome",
                                               value=sequence_not_in_uniprot)
        ]
