from neofox.MHC_predictors.MixMHCpred.mixmhc2pred import MixMHC2pred
from neofox.MHC_predictors.MixMHCpred.mixmhcpred import MixMHCpred
from neofox.MHC_predictors.netmhcpan.combine_netmhcIIpan_pred_multiple_binders import BestAndMultipleBinderMhcII
from neofox.MHC_predictors.netmhcpan.combine_netmhcpan_pred_multiple_binders import BestAndMultipleBinder
from neofox.MHC_predictors.prime import Prime
from neofox.annotation_resources.uniprot.uniprot import Uniprot
from neofox.helpers.blastp_runner import BlastpRunner
from neofox.helpers.runner import Runner
from neofox.model.mhc_parser import MhcParser
from neofox.model.neoantigen import Neoantigen, Patient
from neofox.references.references import DependenciesConfiguration, AvailableAlleles, ReferenceFolder, \
    ORGANISM_HOMO_SAPIENS


class NeoantigenMhcBindingAnnotator:

    def __init__(self, references: ReferenceFolder, configuration: DependenciesConfiguration,
                 uniprot: Uniprot, proteome_blastp_runner: BlastpRunner):
        """class to annotate neoantigens"""
        self.runner = Runner()
        self.configuration = configuration
        self.proteome_db = references.proteome_db
        self.available_alleles = references.get_available_alleles()
        self.organism = references.organism
        self.uniprot = uniprot
        self.proteome_blastp_runner = proteome_blastp_runner

        self.mhc_database = references.get_mhc_database()
        self.mhc_parser = MhcParser.get_mhc_parser(self.mhc_database)

    def get_mhc_binding_annotations(self, neoantigen: Neoantigen, patient: Patient):

        has_mhc1 = patient.mhc1 is not None and len(patient.mhc1) > 0
        has_mhc2 = patient.mhc2 is not None and len(patient.mhc2) > 0

        netmhcpan = None
        netmhc2pan = None
        mixmhcpred = None
        mixmhc2pred = None
        prime = None

        if has_mhc1:
            netmhcpan = self.run_netmhcpan(
                self.runner,
                self.configuration,
                self.available_alleles,
                self.mhc_parser,
                neoantigen,
                patient)
        if has_mhc2:
            netmhc2pan = self._run_netmhc2pan(
                self.runner,
                self.configuration,
                self.available_alleles,
                self.mhc_parser,
                neoantigen,
                patient
            )
        # avoids running MixMHCpred and PRIME for non human organisms
        if self.organism == ORGANISM_HOMO_SAPIENS:
            if self.configuration.mix_mhc2_pred is not None and has_mhc2:
                mixmhc2pred = self._run_mixmhc2pred(
                    self.runner,
                    self.configuration,
                    self.mhc_parser,
                    neoantigen,
                    patient,
                )
            if self.configuration.mix_mhc_pred is not None and has_mhc1:
                mixmhcpred = self._run_mixmhcpred(
                    self.runner,
                    self.configuration,
                    self.mhc_parser,
                    neoantigen,
                    patient,
                )
            if self.configuration.mix_mhc_pred is not None and self.configuration.prime is not None and has_mhc1:
                prime = self._run_prime(
                    self.runner,
                    self.configuration,
                    self.mhc_parser,
                    neoantigen,
                    patient,
                )

        return mixmhc2pred, mixmhcpred, netmhc2pan, netmhcpan, prime

    def run_netmhcpan(
            self,
            runner: Runner,
            configuration: DependenciesConfiguration,
            available_alleles: AvailableAlleles,
            mhc_parser: MhcParser,
            neoantigen: Neoantigen,
            patient: Patient,
    ):
        netmhcpan = BestAndMultipleBinder(runner=runner, configuration=configuration, mhc_parser=mhc_parser,
                                          blastp_runner=self.proteome_blastp_runner)
        netmhcpan.run(
            neoantigen=neoantigen,
            mhc1_alleles_patient=patient.mhc1,
            mhc1_alleles_available=available_alleles.get_available_mhc_i(),
            uniprot=self.uniprot,
        )
        return netmhcpan

    def _run_netmhc2pan(
            self,
            runner: Runner,
            configuration: DependenciesConfiguration,
            available_alleles: AvailableAlleles,
            mhc_parser: MhcParser,
            neoantigen: Neoantigen,
            patient: Patient,
    ):
        netmhc2pan = BestAndMultipleBinderMhcII(
            runner=runner, configuration=configuration, mhc_parser=mhc_parser,
            blastp_runner=self.proteome_blastp_runner)
        netmhc2pan.run(
            neoantigen=neoantigen,
            mhc2_alleles_patient=patient.mhc2,
            mhc2_alleles_available=available_alleles.get_available_mhc_ii(),
            uniprot=self.uniprot
        )
        return netmhc2pan

    def _run_mixmhcpred(
            self,
            runner: Runner,
            configuration: DependenciesConfiguration,
            mhc_parser: MhcParser,
            neoantigen: Neoantigen,
            patient: Patient,
    ):
        mixmhc = MixMHCpred(runner, configuration, mhc_parser)
        mixmhc.run(neoantigen=neoantigen, mhc=patient.mhc1, uniprot=self.uniprot)
        return mixmhc

    def _run_prime(
            self,
            runner: Runner,
            configuration: DependenciesConfiguration,
            mhc_parser: MhcParser,
            neoantigen: Neoantigen,
            patient: Patient,
    ):
        prime = Prime(runner, configuration, mhc_parser)
        prime.run(neoantigen=neoantigen, mhc=patient.mhc1, uniprot=self.uniprot)
        return prime

    def _run_mixmhc2pred(
            self,
            runner: Runner,
            configuration: DependenciesConfiguration,
            mhc_parser: MhcParser,
            neoantigen: Neoantigen,
            patient: Patient,
    ):
        mixmhc2 = MixMHC2pred(runner, configuration, mhc_parser)
        mixmhc2.run(mhc=patient.mhc2, neoantigen=neoantigen, uniprot=self.uniprot)
        return mixmhc2