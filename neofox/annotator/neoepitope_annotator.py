#!/usr/bin/env python
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

from logzero import logger
from datetime import datetime
import neofox
from neofox.MHC_predictors.MixMHCpred.mixmhc2pred import MixMHC2pred
from neofox.MHC_predictors.MixMHCpred.mixmhcpred import MixMHCpred
from neofox.MHC_predictors.prime import Prime
from neofox.annotator.abstract_annotator import AbstractAnnotator
from neofox.annotator.neoantigen_mhc_binding_annotator import NeoantigenMhcBindingAnnotator
from neofox.annotator.neoepitope_mhc_binding_annotator import NeoepitopeMhcBindingAnnotator
from neofox.helpers.epitope_helper import EpitopeHelper
from neofox.MHC_predictors.netmhcpan.combine_netmhcpan_pred_multiple_binders import BestAndMultipleBinder
from neofox.model.factories import AnnotationFactory
from neofox.model.mhc_parser import MhcParser
from neofox.published_features.Tcell_predictor.tcellpredictor_wrapper import TcellPrediction
from neofox.published_features.neoag.neoag_gbm_model import NeoagCalculator
from neofox.published_features.self_similarity.self_similarity import SelfSimilarityCalculator
from neofox.published_features.expression import Expression
from neofox.model.neoantigen import Patient, Neoantigen, Annotations, PredictedEpitope
from neofox.published_features.vaxrank.vaxrank import VaxRank
from neofox.references.references import (
    ReferenceFolder,
    DependenciesConfiguration,
    ORGANISM_HOMO_SAPIENS
)


class NeoepitopeAnnotator(AbstractAnnotator):
    def __init__(self, references: ReferenceFolder, configuration: DependenciesConfiguration,
                 tcell_predictor: TcellPrediction, self_similarity: SelfSimilarityCalculator):
        """class to annotate neoantigens"""

        super().__init__(references, configuration, tcell_predictor, self_similarity)
        self.proteome_db = references.proteome_db
        self.available_alleles = references.get_available_alleles()

        # NOTE: these resources do not read any file thus can be initialised fast
        self.neoag_calculator = NeoagCalculator(runner=self.runner, configuration=configuration)
        self.mhc_database = references.get_mhc_database()
        self.mhc_parser = MhcParser.get_mhc_parser(self.mhc_database)

        self.neoepitope_mhc_binding_annotator = NeoepitopeMhcBindingAnnotator(
            references=references, configuration=configuration, proteome_blastp_runner=self.proteome_blastp_runner,
            uniprot=self.uniprot)

        self.resources_versions = references.get_resources_versions()

    def get_annotated_neoepitope(self, neoepitope: PredictedEpitope) -> PredictedEpitope:
        neoepitope.neofox_annotations = Annotations(
            annotator="NeoFox",
            annotator_version=neofox.VERSION,
            timestamp="{:%Y%m%d%H%M%S%f}".format(datetime.now()),
            resources=self.resources_versions,
            annotations=[]
        )

        # if the WT is not provided it searches for the closest match in the proteome
        if neoepitope.wild_type_peptide is None or neoepitope.wild_type_peptide == '':
            neoepitope.wild_type_peptide = self.proteome_blastp_runner.get_most_similar_wt_epitope(
                neoepitope.mutated_peptide)

        # Runs netmhcpan, netmhc2pan, mixmhcpred and mixmhc2prd in parallel
        annotated_neoepitope = self.neoepitope_mhc_binding_annotator.get_mhc_binding_annotations(neoepitope=neoepitope)

        has_mhc1 = annotated_neoepitope.allele_mhc_i is not None and annotated_neoepitope.allele_mhc_i.name

        if has_mhc1:
            annotated_neoepitope = self.get_additional_annotations_neoepitope_mhci(epitope=annotated_neoepitope)
        else:
            annotated_neoepitope = self.get_additional_annotations_neoepitope_mhcii(epitope=annotated_neoepitope)

        return annotated_neoepitope
