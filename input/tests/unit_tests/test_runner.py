import unittest
from unittest import TestCase

from input.helpers.runner import Runner


class TestRunner(TestCase):

    def setUp(self):
        self.runner = Runner()

    def test_runner(self):
        output, errors = self.runner.run_command(cmd=['python', '-V'])
        self.assertTrue('Python 3.7' in output)
        self.assertTrue(len(errors) == 0)

    def test_runner_failure(self):
        with self.assertRaises(Exception):
            self.runner.run_command(cmd=['nocommandwiththisname'])


if __name__ == "__main__":
    unittest.main()
