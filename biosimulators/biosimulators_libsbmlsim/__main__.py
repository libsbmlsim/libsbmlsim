""" BioSimulators-compliant command-line interface to the `LibSBMLsim <http://fun.bio.keio.ac.jp/software/libsbmlsim/>`_ library.

:Author: Jonathan Karr <karr@mssm.edu>
:Date: 2021-03-27
:Copyright: 2021, BioSimulators Team
:License: MIT
"""

from . import get_simulator_version
from ._version import __version__
from .core import exec_sedml_docs_in_combine_archive
from biosimulators_utils.simulator.cli import build_cli

App = build_cli('libsbmlsim', __version__,
                'LibSBMLsim', get_simulator_version(), 'http://fun.bio.keio.ac.jp/software/libsbmlsim/',
                exec_sedml_docs_in_combine_archive)


def main():
    with App() as app:
        app.run()
