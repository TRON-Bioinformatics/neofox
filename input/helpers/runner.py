import subprocess
import time

from logzero import logger

from input.exceptions import INPuTCommandException


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
