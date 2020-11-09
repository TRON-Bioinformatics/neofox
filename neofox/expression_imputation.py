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


import pysam
from logzero import logger

from neofox.model.neoantigen import Annotation
from neofox.model.wrappers import AnnotationFactory


class ExpressionAnnotator(object):

    def __init__(self, expression_file, tcga_cohort_index_file):
        """
        The expression file is tab separated values file compressed with bgzip and tabix indexed on the gene name and the
        TCGA cohort (encoded as an integer.
        The header of the file is as follows:
        #Gene_Symbol    TCGA_Cohort     Median_Exp      Exp_90p Exp_10p Study_Name
        """
        self.cohort_indices = self._load_tcga_cohort_indices(tcga_cohort_index_file)
        self.expression = pysam.TabixFile(expression_file)

    def _load_tcga_cohort_indices(self, tcga_cohort_index_file):
        data = {}
        with open(tcga_cohort_index_file, "rt") as p:
            for line in p:
                if not line.startswith("TCGA_Cohort"):
                    splitted_line = line.rstrip("\n").split("\t")
                    data[splitted_line[1]] = int(splitted_line[0])
        return data

    def get_gene_expression_annotation(self, gene_name: str, tcga_cohort: str) -> float:
        """
        Returns median gene expression of given gene in TCGA cohort
        :param gene_name: gene name
        :param tcga_cohort: cancer entity, needs to be on of the available TCGA cohorts
        :return: the gene expression
        """
        expression_value = None
        try:
            cohort_index = self.cohort_indices[tcga_cohort]
            expression_value = self._get_gene_expression(gene_name, cohort_index)
        except KeyError:
            logger.error("Tumor type is not available in TCGA data")
        return expression_value

    def _get_gene_expression(self, gene_name: str, cohort_index: int) -> float:
        median_gene_expression = None
        logger.info("Fetching the gene expression at {}:{}".format(gene_name, cohort_index))
        try:
            results = self.expression.fetch(gene_name, cohort_index - 1, cohort_index)
            median_gene_expression = float(next(results).split("\t")[2])
            logger.info("Fetched a gene expression of {}".format(median_gene_expression))
        except (StopIteration, ValueError):
            # pysam triggers these two exceptions under different situations when no data is available
            logger.error("No gene expression entries for {}:{}".format(gene_name, cohort_index))
        except TypeError:
            logger.error("Non float entry coming out of gene expression table for {}:{}".
                         format(gene_name, cohort_index))
        return median_gene_expression
