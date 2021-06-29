""" Tests of the command-line interface

:Author: Jonathan Karr <karr@mssm.edu>
:Date: 2020-10-29
:Copyright: 2020, Center for Reproducible Biomedical Modeling
:License: MIT
"""

from biosimulators_libsbmlsim import __main__
from biosimulators_libsbmlsim import core
from biosimulators_libsbmlsim.data_model import KISAO_ALGORITHMS_MAP
from biosimulators_utils.combine import data_model as combine_data_model
from biosimulators_utils.combine.io import CombineArchiveWriter
from biosimulators_utils.sedml.data_model import (
    SedDocument, Model, ModelLanguage, UniformTimeCourseSimulation, Task, Variable, Symbol,
    Algorithm, AlgorithmParameterChange,
    Report, DataGenerator, DataSet)
from biosimulators_utils.sedml.io import SedmlSimulationWriter
from biosimulators_utils.report import data_model as report_data_model
from biosimulators_utils.report.io import ReportReader
from biosimulators_utils.warnings import BioSimulatorsWarning
from kisao.exceptions import AlgorithmCannotBeSubstitutedException
from kisao.warnings import AlgorithmSubstitutedWarning
from unittest import mock
import biosimulators_libsbmlsim
import copy
import datetime
import dateutil.tz
import os
import numpy.testing
import shutil
import tempfile
import unittest


class CliTestCase(unittest.TestCase):
    FIXTURE = os.path.join(os.path.dirname(__file__), 'fixtures', 'BIOMD0000000075.xml')
    NAMESPACES = {
        'sbml': 'http://www.sbml.org/sbml/level3/version1/core',
    }

    def setUp(self):
        self.dirname = tempfile.mkdtemp()

    def tearDown(self):
        shutil.rmtree(self.dirname)

    def test_version(self):
        with open(os.path.join(os.path.dirname(__file__), '..', '..', 'src', 'libsbmlsim', 'version.h'), 'r') as file:
            for line in file:
                if 'LIBSBMLSIM_DOTTED_VERSION' in line:
                    _, version, _ = line.split('"')
        self.assertEqual(biosimulators_libsbmlsim.__version__,
                         version)

    def test_exec_sed_task(self):
        task = Task(
            model=Model(source=self.FIXTURE, language=ModelLanguage.SBML),
            simulation=UniformTimeCourseSimulation(
                initial_time=0.,
                output_start_time=0.,
                output_end_time=10.,
                number_of_steps=10,
                algorithm=Algorithm(
                    kisao_id='KISAO_0000030',
                )
            )
        )
        variables = [
            Variable(
                id='time',
                symbol=Symbol.time.value,
                task=task,
            ),
            Variable(
                id='PIP2_PHGFP_PM',
                target="/sbml:sbml/sbml:model/sbml:listOfSpecies/sbml:species[@id='PIP2_PHGFP_PM']",
                target_namespaces=self.NAMESPACES,
                task=task,
            ),
        ]
        results, log = core.exec_sed_task(task, variables)
        self.assertEqual(set(results.keys()), set(['time', 'PIP2_PHGFP_PM']))
        numpy.testing.assert_allclose(results['time'], numpy.linspace(0., 10., 11))
        self.assertEqual(results['PIP2_PHGFP_PM'].shape, (11,))
        self.assertFalse(numpy.any(numpy.isnan(results['PIP2_PHGFP_PM'])))

        task2 = copy.deepcopy(task)
        task2.simulation.output_start_time = 5.
        results, log = core.exec_sed_task(task2, variables)
        self.assertEqual(set(results.keys()), set(['time', 'PIP2_PHGFP_PM']))
        numpy.testing.assert_allclose(results['time'], numpy.linspace(5., 10., 11))
        self.assertEqual(results['PIP2_PHGFP_PM'].shape, (11,))
        self.assertFalse(numpy.any(numpy.isnan(results['PIP2_PHGFP_PM'])))

        task2 = copy.deepcopy(task)
        task2.simulation.algorithm.changes.append(AlgorithmParameterChange(
            kisao_id='KISAO_0000483',
            new_value='0.001',
        ))
        results, log = core.exec_sed_task(task2, variables)
        self.assertEqual(set(results.keys()), set(['time', 'PIP2_PHGFP_PM']))
        numpy.testing.assert_allclose(results['time'], numpy.linspace(0., 10., 11))
        self.assertEqual(results['PIP2_PHGFP_PM'].shape, (11,))
        self.assertFalse(numpy.any(numpy.isnan(results['PIP2_PHGFP_PM'])))

        # algorithm substitution
        task2 = copy.deepcopy(task)
        task2.simulation.algorithm.kisao_id = 'KISAO_0000019'
        with mock.patch.dict(os.environ, {'ALGORITHM_SUBSTITUTION_POLICY': 'NONE'}):
            with self.assertRaises(AlgorithmCannotBeSubstitutedException):
                core.exec_sed_task(task2, variables)

        task2 = copy.deepcopy(task)
        task2.simulation.algorithm.kisao_id = 'KISAO_0000019'
        with mock.patch.dict(os.environ, {'ALGORITHM_SUBSTITUTION_POLICY': 'SIMILAR_VARIABLES'}):
            with self.assertWarns(AlgorithmSubstitutedWarning):
                core.exec_sed_task(task2, variables)

        task2 = copy.deepcopy(task)
        task2.simulation.algorithm.changes.append(AlgorithmParameterChange(kisao_id='KISAO_0000488', new_value='1'))
        with mock.patch.dict(os.environ, {'ALGORITHM_SUBSTITUTION_POLICY': 'NONE'}):
            with self.assertRaises(NotImplementedError):
                core.exec_sed_task(task2, variables)

        task2 = copy.deepcopy(task)
        task2.simulation.algorithm.changes.append(AlgorithmParameterChange(kisao_id='KISAO_0000488', new_value='1'))
        with mock.patch.dict(os.environ, {'ALGORITHM_SUBSTITUTION_POLICY': 'SIMILAR_VARIABLES'}):
            with self.assertWarns(BioSimulatorsWarning):
                core.exec_sed_task(task2, variables)

        # all algorithms
        for alg_props in KISAO_ALGORITHMS_MAP.values():
            for order in alg_props['orders']:
                task2 = copy.deepcopy(task)
                task2.simulation.algorithm.kisao_id = alg_props['id']
                task2.simulation.algorithm.changes = []
                if order is not None:
                    task2.simulation.algorithm.changes.append(AlgorithmParameterChange(
                        kisao_id='KISAO_0000594',
                        new_value=str(order),
                    ))
                task2.simulation.algorithm.changes.append(AlgorithmParameterChange(
                    kisao_id='KISAO_0000483',
                    new_value='0.001',
                ))
                results, log = core.exec_sed_task(task2, variables)
                self.assertEqual(set(results.keys()), set(['time', 'PIP2_PHGFP_PM']))
                numpy.testing.assert_allclose(results['time'], numpy.linspace(0., 10., 11))
                self.assertEqual(results['PIP2_PHGFP_PM'].shape, (11,))
                self.assertFalse(numpy.any(numpy.isnan(results['PIP2_PHGFP_PM'])))

        # error handling
        task2 = copy.deepcopy(task)
        task2.simulation.initial_time = 1.
        task2.simulation.output_start_time = 1.
        with self.assertRaisesRegex(NotImplementedError, 'must be zero'):
            core.exec_sed_task(task2, variables)

        task2 = copy.deepcopy(task)
        task2.simulation.algorithm.changes.append(AlgorithmParameterChange(
            kisao_id='KISAO_0000483', new_value='3',
        ))
        with self.assertRaisesRegex(ValueError, 'must specify a positive integer number of steps'):
            core.exec_sed_task(task2, variables)

        task2 = copy.deepcopy(task)
        task2.simulation.output_start_time = 5.01
        task2.simulation.algorithm.changes.append(AlgorithmParameterChange(
            kisao_id='KISAO_0000483', new_value='0.1',
        ))
        with self.assertRaisesRegex(ValueError, 'must specify a positive integer number of intervals'):
            core.exec_sed_task(task2, variables)

        task2 = copy.deepcopy(task)
        task2.model.source = 'undefined'
        with mock.patch('biosimulators_utils.sedml.validation.validate_variable_xpaths', return_value={}):
            with mock.patch('biosimulators_utils.sedml.validation.validate_model', return_value=([], [])):
                with self.assertRaisesRegex(ValueError, 'File Not Found'):
                    core.exec_sed_task(task2, variables)

        variables2 = copy.deepcopy(variables)
        variables2[0].symbol = 'time'
        with self.assertRaisesRegex(NotImplementedError, 'unsupported symbols'):
            core.exec_sed_task(task, variables2)

        variables2 = copy.deepcopy(variables)
        variables2[1].target = '/sbml:sbml/sbml:model'
        with self.assertRaisesRegex(NotImplementedError, 'unsupported targets'):
            core.exec_sed_task(task, variables2)

    def test_exec_sedml_docs_in_combine_archive(self):
        archive_dirname = os.path.join(self.dirname, 'archive')
        os.mkdir(archive_dirname)

        task = Task(
            id='task',
            model=Model(id='model', source=self.FIXTURE, language=ModelLanguage.SBML),
            simulation=UniformTimeCourseSimulation(
                id='simulation',
                initial_time=0.,
                output_start_time=0.,
                output_end_time=10.,
                number_of_steps=10,
                algorithm=Algorithm(
                    kisao_id='KISAO_0000030',
                )
            )
        )
        variables = [
            Variable(
                id='time',
                symbol=Symbol.time.value,
                task=task,
            ),
            Variable(
                id='PIP2_PHGFP_PM',
                target="/sbml:sbml/sbml:model/sbml:listOfSpecies/sbml:species[@id='PIP2_PHGFP_PM']",
                target_namespaces=self.NAMESPACES,
                task=task,
            ),
        ]

        doc = SedDocument(
            models=[task.model],
            simulations=[task.simulation],
            tasks=[task],
        )
        report = Report(id='report')
        doc.outputs.append(report)
        for variable in variables:
            data_gen = DataGenerator(
                id='data_generator_' + variable.id,
                variables=[variable],
                math=variable.id,
            )
            doc.data_generators.append(data_gen)

            report.data_sets.append(DataSet(
                id='data_set_' + variable.id,
                label=variable.id,
                data_generator=data_gen,
            ))

        model_filename = os.path.join(archive_dirname, 'model.xml')
        shutil.copyfile(self.FIXTURE, model_filename)

        sim_filename = os.path.join(archive_dirname, 'sim.sedml')
        SedmlSimulationWriter().run(doc, sim_filename)

        updated = datetime.datetime(2020, 1, 2, 1, 2, 3, tzinfo=dateutil.tz.tzutc())
        archive = combine_data_model.CombineArchive(
            contents=[
                combine_data_model.CombineArchiveContent(
                    'model.xml', combine_data_model.CombineArchiveContentFormat.SBML.value, updated=updated),
                combine_data_model.CombineArchiveContent(
                    'sim.sedml', combine_data_model.CombineArchiveContentFormat.SED_ML.value, updated=updated),
            ],
            updated=updated,
        )
        archive_filename = os.path.join(self.dirname, 'archive.omex')
        CombineArchiveWriter().run(archive, archive_dirname, archive_filename)

        out_dir = os.path.join(self.dirname, 'results')
        core.exec_sedml_docs_in_combine_archive(
            archive_filename, out_dir,
            report_formats=[report_data_model.ReportFormat.h5])

        results = ReportReader().run(report, out_dir, 'sim.sedml/report', format=report_data_model.ReportFormat.h5)

        numpy.testing.assert_allclose(results['data_set_time'], numpy.linspace(0., 10., 11))
        self.assertEqual(results['data_set_PIP2_PHGFP_PM'].shape, (11,))
        self.assertFalse(numpy.any(numpy.isnan(results['data_set_PIP2_PHGFP_PM'])))

    def test_raw_cli(self):
        with mock.patch('sys.argv', ['', '--help']):
            with self.assertRaises(SystemExit) as context:
                __main__.main()
                self.assertRegex(context.Exception, 'usage: ')
