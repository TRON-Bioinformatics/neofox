from math import ceil, floor

from Bio import SeqIO
from Bio.Align import substitution_matrices
from Bio.Alphabet.IUPAC import ExtendedIUPACProtein


class PyHex:

    def __init__(self, iedb_fasta, magic_number=4):
        self.iedb_sequences = self._read_fasta(iedb_fasta)
        self.magic_number = magic_number
        self.blosum = substitution_matrices.load("BLOSUM62")

    @staticmethod
    def _read_fasta(fasta_file):
        sequences = []
        # read fasta
        with open(fasta_file, "r") as handle:
            for record in SeqIO.parse(handle, "fasta"):
                # include only records that do not contain non-standard amino acids
                if not any([aa not in ExtendedIUPACProtein.letters for aa in record.seq]):
                    sequences.append(record)
        return sequences

    def _align(self, sequence, mutated_sequence):
        weights = self._get_sequence_weights(mutated_sequence)
        score = sum([self.blosum[q, t] * w for q, t, w in zip(sequence, mutated_sequence, weights)])
        return score

    def _get_sequence_weights(self, mutated_sequence):
        length_mutated_sequence = len(mutated_sequence)
        mid_score = ceil(length_mutated_sequence / 2) * self.magic_number
        weights = list(range(1, mid_score, self.magic_number))
        weights.extend(reversed(weights[0:floor(length_mutated_sequence / 2)]))

        top_floor = floor(length_mutated_sequence / 3)
        weights[0:top_floor] = list(range(1, top_floor + 1))
        tail = length_mutated_sequence - top_floor
        weights[tail:length_mutated_sequence] = list(reversed(range(1, top_floor + 1)))

        return weights

    def run(self, mutated_sequence):
        # excludes sequences that have different length than the mutated sequence
        sequences = [s for s in self.iedb_sequences if len(s.seq) == len(mutated_sequence)]
        # align each of the sequences
        alignment_scores = [self._align(s.seq, mutated_sequence) for s in sequences]
        # gets the best score of all the alignments
        best_score = max(alignment_scores)
        return best_score
