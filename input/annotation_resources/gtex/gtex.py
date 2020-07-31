import os
import pandas as pd
from logzero import logger


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

    def _get_metric(self, gene, tissue, metric):
        # TODO: avoid the stupid NA
        try:
            metric = self.gtex.at[(gene, tissue, metric), 'value']
        except KeyError:
            metric = 'NA'
        return metric