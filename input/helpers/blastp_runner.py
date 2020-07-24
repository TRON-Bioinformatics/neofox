from input.helpers import intermediate_files
from input.predictors.neoantigen_fitness.aligner import Aligner


class BlastpRunner(object):

    def __init__(self, runner, configuration):
        """
        :type runner: input.helpers.runner.Runner
        :type configuration: input.references.DependenciesConfiguration
        """
        self.runner = runner
        self.configuration = configuration

    def run_blastp(self, fasta_file, database):
        '''
        This function runs BLASTP on a given database
        '''
        outfile = intermediate_files.create_temp_file(prefix="tmp_blastp_", suffix=".xml")
        self.runner.run_command(cmd=[
            self.configuration.blastp,
            "-gapopen", "11",
            "-gapextend", "1",
            "-outfmt", "5",
            "-query", fasta_file,
            "-out", outfile,
            "-db", database,
            "-evalue", "100000000"])
        return outfile

    def parse_blastp_output(self, blastp_output_file, **kwargs):
        aligner = Aligner()
        # set a to 32 for dissimilarity
        aligner.readAllBlastAlignments(blastp_output_file)
        aligner.computeR(**kwargs)
        return aligner.Ri.get(1, 0)  # NOTE: returns 0 when not present