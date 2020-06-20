import re

from logzero import logger


class ProveanAnnotator(object):

    def __init__(self, provean_file, header_epitopes, epitopes):
        """
        Loads provean scores as dictionary, but only for ucsc ids that are in epitope list
        """
        epitope_ids = self._load_ucsc_ids_epitopes(header_epitopes=header_epitopes, epitopes=epitopes)
        logger.info("Starting load of PROVEAN matrix" + provean_file)
        self.header_provean, self.provean_matrix = self._load_provean_matrix(epitope_ids, provean_file)
        logger.info("PROVEAN matrix loaded")

    def _load_ucsc_ids_epitopes(self, header_epitopes, epitopes):
        """
        Returns set with ucsc ids of epitopes.
        """
        col_ucsc = header_epitopes.index("UCSC_transcript")
        col_pos = header_epitopes.index("substitution")
        return set([self.build_ucsc_id_plus_position(ucsc_id=e[col_ucsc], substitution=e[col_pos]) for e in epitopes])

    def _load_provean_matrix(self, epitope_ids, provean_file):
        provean_matrix = {}
        with open(provean_file) as f:
            header = next(f).rstrip().split(";")  # stores header
            for line in f:
                parts = line.rstrip().split(";")
                ucsc_id_pos = parts[-1]
                if ucsc_id_pos in epitope_ids:
                    provean_matrix[ucsc_id_pos] = parts
        return header, provean_matrix

    def get_provean_annotation(self, mutated_aminoacid, ucsc_id_position):
        """
        This function maps Provean score on given position and for specific SNV onto epitope data set
        (which is in form of tuple --> header + dict of ucsc_pos_id: df row)
        """
        try:
            return self.provean_matrix[ucsc_id_position][self.header_provean.index(mutated_aminoacid)]
        except (ValueError, KeyError) as e:
            return "NA"

    @staticmethod
    def build_ucsc_id_plus_position(substitution, ucsc_id):
        ucsc_epi = re.sub(r'.\d+$', '', ucsc_id)
        position_match = re.match(r'[A-Z](\d+)[A-Z]', substitution)
        pos_prot = position_match.group(1) if position_match else "Del"
        return "{}_{}".format(ucsc_epi, pos_prot)
