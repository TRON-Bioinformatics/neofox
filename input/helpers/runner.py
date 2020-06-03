import subprocess
import time

from logzero import logger

from input.exceptions import INPuTCommandException
from input.helpers import intermediate_files
from input.neoantigen_fitness.aligner import Aligner


class Runner(object):

    def run_command(self, cmd, **kwargs):
        logger.info("Starting command: {}".format(" ".join(cmd)))
        start = time.time()
        process = subprocess.Popen(self._preprocess_command(cmd), stderr=subprocess.PIPE, stdout=subprocess.PIPE, **kwargs)
        output, errors = process.communicate()
        return_code = process.returncode
        end = time.time()
        logger.info("Elapsed time {} seconds".format(int(end - start)))
        if return_code == 0:
            logger.info("Finished command correctly!")
            logger.info(self._decode(output))
        else:
            logger.error("Finished command with return code {}".format(return_code))
            logger.error(self._decode(output))
            logger.error(self._decode(errors))
            raise INPuTCommandException("Error running command '{}'".format(" ".join(cmd)))
        return self._decode(output), self._decode(errors)

    @staticmethod
    def _preprocess_command(cmd):
        """
        This makes sure that any parameter containing white spaces is passed appropriately
        """
        return " ".join(cmd).split(" ")

    def _decode(self, data):
        return data.decode('utf8')


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
