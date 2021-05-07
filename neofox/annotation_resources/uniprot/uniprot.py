#
# Copyright (c) 2020-2030 Translational Oncology at the Medical Center of the Johannes Gutenberg-University Mainz gGmbH.
#
# This file is part of Neofox
# (see https://github.com/tron-bioinformatics/neofox).
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program. If not, see <http://www.gnu.org/licenses/>.#
import pickle
from typing import List
from logzero import logger
from neofox.model.neoantigen import Annotation
from neofox.model.wrappers import AnnotationFactory


class Uniprot(object):
    """
    Loads the whole Uniprot fasta in a single string to check for a exact match of an aminoacid sequence.
    Even though this is not the most elegant, this is the fastest by 1 order of magnitude in both
    data loading and checking for an exact match
    """

    def __init__(self, proteome):
        logger.info("Loading Uniprot...")
        self.uniprot = self._load_proteome(proteome)
        logger.info("Loaded Uniprot.")

    @staticmethod
    def _load_proteome(proteome) -> str:
        with open(proteome, 'rb') as f:
            uniprot_unpickled = pickle.load(f)
        return uniprot_unpickled

    def is_sequence_not_in_uniprot(self, sequence) -> bool:
        return self.uniprot.find(sequence) < 0

    def get_annotations(self, sequence_not_in_uniprot: bool) -> List[Annotation]:
        return [
            AnnotationFactory.build_annotation(
                name="mutation_not_found_in_proteome", value=sequence_not_in_uniprot
            )
        ]
