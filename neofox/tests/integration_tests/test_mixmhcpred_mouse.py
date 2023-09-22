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
from unittest import TestCase
from logzero import logger

from neofox.model.factories import MhcFactory
from neofox.model.mhc_parser import MhcParser
from neofox.model.neoantigen import Mhc2Name
from neofox.helpers.epitope_helper import EpitopeHelper
import neofox.tests.integration_tests.integration_test_tools as integration_test_tools
from neofox.MHC_predictors.MixMHCpred.mixmhc2pred import MixMHC2pred
from neofox.helpers.runner import Runner
from neofox.annotation_resources.uniprot.uniprot import Uniprot
from neofox.tests.tools import get_neoantigen
from neofox.references.references import ReferenceFolder, DependenciesConfiguration, ORGANISM_HOMO_SAPIENS, \
    ORGANISM_MUS_MUSCULUS


class TestMixMHCPredMouse(TestCase):
    def setUp(self):
        self.references, self.configuration = integration_test_tools.load_references(organism=ORGANISM_MUS_MUSCULUS)
        self.runner = Runner()
        mhc_parser = MhcParser.get_mhc_parser(self.references.get_mhc_database())
        #self.mixmhcpred = MixMHCpred(
        #    runner=self.runner, configuration=self.configuration, mhc_parser=mhc_parser
        #)
        self.mixmhc2pred = MixMHC2pred(
            runner=self.runner, configuration=self.configuration, mhc_parser=mhc_parser,
            references=self.references
        )
        self.hla_database = self.references.get_mhc_database()
        #self.test_mhc_one = integration_test_tools.get_hla_one_test(self.hla_database)
        self.test_mhc_two = integration_test_tools.get_h2_two_test(self.hla_database)
        self.test_mhc_two_b = integration_test_tools.get_h2_two_test_b(self.hla_database)
        self.uniprot = Uniprot(self.references.uniprot_pickle)

    def test_mixmhcpred2_antigen_iedb_b_haplotype(self):
        # Test mixmhc2pred with H2Ab allele (C57BL/6 setting)
        # this is an antigen from IEDB of length 27
        neoantigen = get_neoantigen(
            mutated_xmer="RQHSIKEGLQFIQPPLSYPGTQEQYAV",
            wild_type_xmer= "RQHSIKEGLQFIQSPLSYPGTQEQYAV")
        self.mixmhc2pred.run(
            neoantigen=neoantigen, mhc=self.test_mhc_two_b,
            uniprot=self.uniprot
        )

        best_result = EpitopeHelper.select_best_by_rank(predictions=self.mixmhc2pred.results)

        self.assertEquals("QPPLSYPGTQEQYAV", best_result.mutated_peptide)
        self.assertEquals(9.43, best_result.rank_mutated)
        self.assertEquals("H2Ab", best_result.isoform_mhc_i_i.name)

    def test_mixmhcpred2_antigen_iedb(self):
        # Test mixmhc2pred with H2Ad and H2Ed allele (BALB/c setting)
        # this is an antigen from IEDB of length 27
        neoantigen = get_neoantigen(
            mutated_xmer="RQHSIKEGLQFIQPPLSYPGTQEQYAV",
            wild_type_xmer= "RQHSIKEGLQFIQSPLSYPGTQEQYAV")
        self.mixmhc2pred.run(
            neoantigen=neoantigen, mhc=self.test_mhc_two,
            uniprot=self.uniprot
        )

        best_result = EpitopeHelper.select_best_by_rank(predictions=self.mixmhc2pred.results)

        self.assertEquals("KEGLQFIQPPLSYPG", best_result.mutated_peptide)
        self.assertEquals(11.9, best_result.rank_mutated)
        self.assertEquals("H2Ad", best_result.isoform_mhc_i_i.name)

    def test_mixmhcpred2_no_mutation(self):
        neoantigen = get_neoantigen(
            mutated_xmer="RQHSIKEGLQFIQSPLSYPGTQEQYAV",
            wild_type_xmer= "RQHSIKEGLQFIQSPLSYPGTQEQYAV")
        self.mixmhc2pred.run(
            neoantigen=neoantigen, mhc=self.test_mhc_two,
            uniprot=self.uniprot
        )

        best_result = EpitopeHelper.select_best_by_rank(predictions=self.mixmhc2pred.results)

        self.assertIsNone(best_result.mutated_peptide)
        self.assertIsNone(best_result.rank_mutated)
        self.assertIsNone(best_result.isoform_mhc_i_i.name)

    def test_mixmhc2pred_allele(self):
        neoantigen = get_neoantigen(mutated_xmer="RQHSIKEGLQFIQPPLSYPGTQEQYAV", wild_type_xmer="RQHSIKEGLQFIQSPLSYPGTQEQYAV")
        # this is a MHC II genotype which results in no available alleles for MixMHC2pred
        MHC_TWO_NEW = MhcFactory.build_mhc2_alleles(
            [
                "H2Ab",
                "H2Ad",
                "H2Ed"
                # this mouse allele is supported by MixMHC2pred but does not exist in H2 database
                #"H2Anb1"
            ],
            self.hla_database
        )
        alleles = self.mixmhc2pred.transform_h2_alleles_for_prediction(MHC_TWO_NEW)
        logger.info(alleles)
        self.assertListEqual(alleles, ['H2_Aa_b__H2_Ab_b', 'H2_Aa_d__H2_Ab_d', 'H2_Ea_d__H2_Eb_d'])
        self.mixmhc2pred.run(neoantigen=neoantigen, mhc=MHC_TWO_NEW, uniprot=self.uniprot)
        best_result = EpitopeHelper.select_best_by_rank(predictions=self.mixmhc2pred.results)
        self.assertIsNotNone(best_result.mutated_peptide)
        self.assertIsNotNone(best_result.rank_mutated)
        self.assertIsNotNone(best_result.isoform_mhc_i_i.name)
