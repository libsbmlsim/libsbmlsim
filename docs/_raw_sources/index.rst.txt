BioSimulators-LibSBMLSim documentation
======================================

BioSimulators-LibSBMLSim provides a `BioSimulators <https://biosimulators.org>`_-compliant command-line interface to the `LibSBMLSim <https://fun.bio.keio.ac.jp/software/libsbmlsim/>`_ simulation tool. A Docker image for this package is also available.

This command-line interface and Docker image enable users to use LibSBMLSim to execute `COMBINE/OMEX archives <https://combinearchive.org/>`_ that describe one or more simulation experiments (in `SED-ML format <https://sed-ml.org>`_) of one or more kinetic models in SBML format.

A list of the algorithms and algorithm parameters supported by LibSBMLSim is available at `BioSimulators <https://biosimulators.org/simulators/libsbmlsim>`_.

A simple web application and web service for using LibSBMLSim to execute COMBINE/OMEX archives is also available at `runBioSimulations <https://run.biosimulations.org>`_.

Contents
--------

.. toctree::
   :maxdepth: 2

   installation.rst
   tutorial.rst
   API documentation <source/biosimulators_libsbmlsim.rst>
   about.rst
   genindex.rst
