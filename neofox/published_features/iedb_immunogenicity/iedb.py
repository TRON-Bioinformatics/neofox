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
from typing import List, Union
from logzero import logger
from neofox.model.neoantigen import Annotation, MhcAllele, PredictedEpitope
from neofox.model.factories import AnnotationFactory

immunoweight = [0.00, 0.00, 0.10, 0.31, 0.30, 0.29, 0.26, 0.18, 0.00]

immunoscale = {
    "A": 0.127,
    "C": -0.175,
    "D": 0.072,
    "E": 0.325,
    "F": 0.380,
    "G": 0.110,
    "H": 0.105,
    "I": 0.432,
    "K": -0.700,
    "L": -0.036,
    "M": -0.570,
    "N": -0.021,
    "P": -0.036,
    "Q": -0.376,
    "R": 0.168,
    "S": -0.537,
    "T": 0.126,
    "V": 0.134,
    "W": 0.719,
    "Y": -0.012,
}

allele_dict = {
    "H-2-Db": "2,5,9",
    "H-2-Dd": "2,3,5",
    "H-2-Kb": "2,3,9",
    "H-2-Kd": "2,5,9",
    "H-2-Kk": "2,8,9",
    "H-2-Ld": "2,5,9",
    "HLA-A*01:01": "2,3,9",
    "HLA-A*02:01": "1,2,9",
    "HLA-A*02:02": "1,2,9",
    "HLA-A*02:03": "1,2,9",
    "HLA-A*02:06": "1,2,9",
    "HLA-A*02:11": "1,2,9",
    "HLA-A*03:01": "1,2,9",
    "HLA-A*11:01": "1,2,9",
    "HLA-A*23:01": "2,7,9",
    "HLA-A*24:02": "2,7,9",
    "HLA-A*26:01": "1,2,9",
    "HLA-A*29:02": "2,7,9",
    "HLA-A*30:01": "1,3,9",
    "HLA-A*30:02": "2,7,9",
    "HLA-A*31:01": "1,2,9",
    "HLA-A*32:01": "1,2,9",
    "HLA-A*33:01": "1,2,9",
    "HLA-A*68:01": "1,2,9",
    "HLA-A*68:02": "1,2,9",
    "HLA-A*69:01": "1,2,9",
    "HLA-B*07:02": "1,2,9",
    "HLA-B*08:01": "2,5,9",
    "HLA-B*15:01": "1,2,9",
    "HLA-B*15:02": "1,2,9",
    "HLA-B*18:01": "1,2,9",
    "HLA-B*27:05": "2,3,9",
    "HLA-B*35:01": "1,2,9",
    "HLA-B*39:01": "1,2,9",
    "HLA-B*40:01": "1,2,9",
    "HLA-B*40:02": "1,2,9",
    "HLA-B*44:02": "2,3,9",
    "HLA-B*44:03": "2,3,9",
    "HLA-B*45:01": "1,2,9",
    "HLA-B*46:01": "1,2,9",
    "HLA-B*51:01": "1,2,9",
    "HLA-B*53:01": "1,2,9",
    "HLA-B*54:01": "1,2,9",
    "HLA-B*57:01": "1,2,9",
    "HLA-B*58:01": "1,2,9",
}


class IEDBimmunogenicity:

    def predict_immunogenicity(self, pep, allele):

        custom_mask = allele_dict.get(allele, False)
        peptide = pep.upper()
        peplen = len(peptide)

        cterm = peplen - 1
        score = 0
        count = 0

        if not custom_mask:
            mask_num = [0, 1, cterm]
        elif custom_mask:
            mask_str = custom_mask.split(",")
            mask_num = list(map(int, mask_str))
            mask_num = [x - 1 for x in mask_num]

        if peplen > 9:
            pepweight = immunoweight[:5] + ((peplen - 9) * [0.30]) + immunoweight[5:]
        else:
            pepweight = immunoweight
        try:
            for pos in peptide:
                if pos not in list(immunoscale.keys()):
                    logger.debug("{} {} {}".format(pos, pep, allele))
                    raise KeyError()
                elif count not in mask_num:
                    score += pepweight[count] * immunoscale[pos]
                    count += 1
                else:
                    count += 1
        except Exception as ex:
            logger.exception(ex)
        return score

    def calculate_iedb_immunogenicity(
        self, peptide, mhc_allele: Union[MhcAllele, None], mhc_score
    ):
        """This function determines the IEDB immunogenicity score"""
        score = None
        try:
            if peptide != "-":
                score = self.predict_immunogenicity(peptide, mhc_allele.name if mhc_allele else None)
                logger.info(score)
        except (ValueError, AttributeError):
            pass
        return score

    def get_annotations(
            self, mutated_peptide_mhci: PredictedEpitope, mutated_peptide_mhcii: PredictedEpitope) -> List[Annotation]:
        """
        returns IEDB immunogenicity for best predicted MHC I (based on affinity) and MHC II (based on affinity) epitopes
        """
        iedb = None
        iedb_mhcii = None
        if mutated_peptide_mhci and mutated_peptide_mhci.mutated_peptide:
            iedb = self.calculate_iedb_immunogenicity(
                        peptide=mutated_peptide_mhci.mutated_peptide,
                        mhc_allele=mutated_peptide_mhci.allele_mhc_i,
                        mhc_score=mutated_peptide_mhci.affinity_mutated,
                    )
        if mutated_peptide_mhcii and mutated_peptide_mhcii.mutated_peptide:
            iedb_mhcii = self.calculate_iedb_immunogenicity(
                peptide=mutated_peptide_mhcii.mutated_peptide,
                mhc_allele=None,
                mhc_score=mutated_peptide_mhcii.affinity_mutated,
            )
        annotations = [
            AnnotationFactory.build_annotation(
                value=iedb,
                name="IEDB_Immunogenicity_MHCI",
            ),
            AnnotationFactory.build_annotation(
                name="IEDB_Immunogenicity_MHCII",
                value=iedb_mhcii
            )
        ]
        return annotations

    def get_annotations_epitope_mhcii(self, epitope: PredictedEpitope) -> List[Annotation]:
        return [
            AnnotationFactory.build_annotation(
                value=self.calculate_iedb_immunogenicity(
                    peptide=epitope.mutated_peptide, mhc_allele=None, mhc_score=epitope.affinity_mutated),
                name='IEDB_Immunogenicity')
        ]

    def get_annotations_epitope_mhci(self, epitope: PredictedEpitope) -> List[Annotation]:
        return [
            AnnotationFactory.build_annotation(
                value=self.calculate_iedb_immunogenicity(
                    peptide=epitope.mutated_peptide, mhc_allele=epitope.allele_mhc_i,
                    mhc_score=epitope.affinity_mutated),
                name='IEDB_Immunogenicity')
        ]
