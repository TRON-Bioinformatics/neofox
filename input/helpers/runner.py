import subprocess

from logzero import logger


class Runner(object):

    def run_command(self, cmd, **kwargs):
        logger.info("Starting command: {}".format(" ".join(cmd)))
        process = subprocess.Popen(cmd, stderr=subprocess.PIPE, stdout=subprocess.PIPE, **kwargs)
        output, errors = process.communicate()
        return_code = process.returncode
        if return_code == 0:
            logger.info("Finished command correctly!")
        else:
            logger.error("Finished command with return code {}".format(return_code))
        return self._decode(output), self._decode(errors)

    def _decode(self, data):
        return data.decode('utf8')
