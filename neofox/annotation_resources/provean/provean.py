import re
import pysam
import gzip
from logzero import logger

from neofox.model.neoantigen import Annotation
from neofox.model.wrappers import AnnotationFactory


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

    def get_provean_annotation(self, protein_id: str, position: int, mutated_aminoacid: str) -> Annotation:
        """
        Returns the PROVEAN score annotation of particular mutation in a protein
        :param protein_id: ucsc protein id
        :param position: the position in the protein
        :param mutated_aminoacid: the mutated aminoacid
        :return: the provean score
        """
        score = self._get_provean_score(mutated_aminoacid, position, protein_id)
        return AnnotationFactory.build_annotation(name="PROVEAN_score", value=score)

    def _get_provean_score(self, mutated_aminoacid:str, position:int, protein_id:str) -> float:
        provean_score = None
        logger.info("Fetching the PROVEAN score at {}:{}:{}".format(protein_id, position, mutated_aminoacid))
        try:
            results = self.provean.fetch(protein_id, position - 1, position)
            provean_entry = next(results).split("\t")
            provean_score = float(provean_entry[self.aminoacid_indexes.get(mutated_aminoacid)])
            logger.info("Fetched a PROVEAN score of {}".format(provean_score))
        except (StopIteration, ValueError):
            # pysam triggers these two exceptions under different situations when no data is available
            logger.error("No PROVEAN entries for {}:{}".format(protein_id, position))
        except IndexError:
            logger.error("Bad aminoacid trying to read PROVEAN {} with index {}".format(
                mutated_aminoacid, self.aminoacid_indexes[mutated_aminoacid]))
        except TypeError:
            logger.error("Non float entry coming out of PROVEAN for {}:{}".format(protein_id, position))
        return provean_score
