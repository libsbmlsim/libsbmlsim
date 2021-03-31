""" Tests of the command-line interface

:Author: Jonathan Karr <karr@mssm.edu>
:Date: 2020-10-29
:Copyright: 2020, Center for Reproducible Biomedical Modeling
:License: MIT
"""

from biosimulators_libsbmlsim import __main__
# from biosimulators_libsbmlsim import core
from unittest import mock
import shutil
import tempfile
import unittest


class CliTestCase(unittest.TestCase):
    DOCKER_IMAGE = 'ghcr.io/libsbmlsim/libsbmlsim:latest'
    NAMESPACES = {
        'sbml': 'http://www.sbml.org/sbml/level2/version4',
    }

    def setUp(self):
        self.dirname = tempfile.mkdtemp()

    def tearDown(self):
        shutil.rmtree(self.dirname)

    def test_raw_cli(self):
        with mock.patch('sys.argv', ['', '--help']):
            with self.assertRaises(SystemExit) as context:
                __main__.main()
                self.assertRegex(context.Exception, 'usage: ')
