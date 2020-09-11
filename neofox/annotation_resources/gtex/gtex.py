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
import os
from typing import List

import pandas as pd
from logzero import logger

from neofox.model.neoantigen import Annotation
from neofox.model.wrappers import AnnotationFactory

GTEX_FILE_NAME = "gtex_combined.csv.gz"


class GTEx(object):

    def __init__(self):

        gtex_file = os.path.join(os.path.abspath(os.path.dirname(__file__)), GTEX_FILE_NAME)
        logger.info("Loading GTEx...")
        self.gtex = self._load_gtex(gtex_file)
        logger.info("GTEx loaded. {} genes, {} tissues".format(
            len(self.gtex.index.get_level_values('gene').unique()),
            len(self.gtex.index.get_level_values('tissue').unique())))

    def _load_gtex(self, gtex_file):
        """
        Loads GTEx resource into a pandas dataframe with columns: gene, tissue, metric, value
        :rtype: pd.DataFrame
        """
        gtex = pd.read_csv(gtex_file, sep=';', decimal=',', compression='gzip')
        del gtex['Unnamed: 0']
        gtex_stacked = gtex.set_index('gene').stack()
        gtex_stacked.index.names = ['gene', 'name']
        gtex_stacked_df = gtex_stacked.to_frame().reset_index()
        gtex_stacked_df['value'] = gtex_stacked_df[0]
        del gtex_stacked_df[0]
        gtex_stacked_df['tissue'] = gtex_stacked_df['name'].transform(lambda x: x.split('.')[0])
        gtex_stacked_df['metric'] = gtex_stacked_df['name'].transform(lambda x: x.split('.')[1])
        del gtex_stacked_df['name']
        gtex_stacked_df.set_index(['gene', 'tissue', 'metric'], inplace=True)
        return gtex_stacked_df

    def get_metrics(self, gene, tissue):
        """
        Gets metrics from GTEx for a given gene in a given tissue
        :param gene: the gene symbol
        :type gene: str
        :param tissue: the tissue
        :type tissue: str
        :return: mean expression, sum expression, standard deviation expression
        :rtype: float, float, float
        """
        mean_expression = self._get_metric(gene, tissue, 'mean')
        sum_expression = self._get_metric(gene, tissue, 'sum')
        sd_expression = self._get_metric(gene, tissue, 'sd')
        return mean_expression, sum_expression, sd_expression

    def get_annotations(self, mean: float, sd: float, sum: float) -> List[Annotation]:
        return [
            AnnotationFactory.build_annotation(name="Reference_expression_mean", value=mean),
            AnnotationFactory.build_annotation(name="Reference_expression_standard_deviation", value=sd),
            AnnotationFactory.build_annotation(name="Reference_expression_mean_sum", value=sum)
        ]

    def _get_metric(self, gene, tissue, metric):
        try:
            metric = self.gtex.at[(gene, tissue, metric), 'value']
        except KeyError:
            metric = None
        return metric