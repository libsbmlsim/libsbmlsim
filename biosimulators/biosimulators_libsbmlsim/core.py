""" Methods for executing SED tasks in COMBINE archives and saving their outputs

:Author: Jonathan Karr <karr@mssm.edu>
:Date: 2021-03-27
:Copyright: 2021, Center for Reproducible Biomedical Modeling
:License: MIT
"""

from .data_model import KISAO_ALGORITHMS_MAP, get_integrator
from biosimulators_utils.combine.exec import exec_sedml_docs_in_archive
from biosimulators_utils.config import get_config
from biosimulators_utils.log.data_model import CombineArchiveLog, TaskLog  # noqa: F401
from biosimulators_utils.viz.data_model import VizFormat  # noqa: F401
from biosimulators_utils.report.data_model import ReportFormat, VariableResults, SedDocumentResults  # noqa: F401
from biosimulators_utils.sedml.data_model import (Task, ModelLanguage, UniformTimeCourseSimulation,  # noqa: F401
                                                  Algorithm, Variable, Symbol)
from biosimulators_utils.sedml import validation
from biosimulators_utils.sedml.exec import exec_sed_doc
from biosimulators_utils.simulator.utils import get_algorithm_substitution_policy
from biosimulators_utils.utils.core import raise_errors_warnings
from biosimulators_utils.warnings import warn, BioSimulatorsWarning
from kisao.data_model import AlgorithmSubstitutionPolicy, ALGORITHM_SUBSTITUTION_POLICY_LEVELS
from kisao.utils import get_preferred_substitute_algorithm_by_ids
import functools
import libsbmlsim
import os
import pandas
import tempfile


__all__ = ['exec_sedml_docs_in_combine_archive', 'exec_sed_task']


def exec_sedml_docs_in_combine_archive(archive_filename, out_dir,
                                       return_results=False,
                                       report_formats=None, plot_formats=None,
                                       bundle_outputs=None, keep_individual_outputs=None,
                                       raise_exceptions=True):
    """ Execute the SED tasks defined in a COMBINE/OMEX archive and save the outputs

    Args:
        archive_filename (:obj:`str`): path to COMBINE/OMEX archive
        out_dir (:obj:`str`): path to store the outputs of the archive

            * CSV: directory in which to save outputs to files
              ``{ out_dir }/{ relative-path-to-SED-ML-file-within-archive }/{ report.id }.csv``
            * HDF5: directory in which to save a single HDF5 file (``{ out_dir }/reports.h5``),
              with reports at keys ``{ relative-path-to-SED-ML-file-within-archive }/{ report.id }`` within the HDF5 file

        return_results (:obj:`bool`, optional): whether to return the result of each output of each SED-ML file
        report_formats (:obj:`list` of :obj:`ReportFormat`, optional): report format (e.g., csv or h5)
        plot_formats (:obj:`list` of :obj:`VizFormat`, optional): report format (e.g., pdf)
        bundle_outputs (:obj:`bool`, optional): if :obj:`True`, bundle outputs into archives for reports and plots
        keep_individual_outputs (:obj:`bool`, optional): if :obj:`True`, keep individual output files
        raise_exceptions (:obj:`bool`, optional): whether to raise exceptions

    Returns:
        :obj:`tuple`:

            * :obj:`SedDocumentResults`: results
            * :obj:`CombineArchiveLog`: log
    """
    sed_doc_executer = functools.partial(exec_sed_doc, exec_sed_task)
    return exec_sedml_docs_in_archive(sed_doc_executer, archive_filename, out_dir,
                                      apply_xml_model_changes=True,
                                      return_results=return_results,
                                      report_formats=report_formats,
                                      plot_formats=plot_formats,
                                      bundle_outputs=bundle_outputs,
                                      keep_individual_outputs=keep_individual_outputs,
                                      raise_exceptions=raise_exceptions)


def exec_sed_task(task, variables, log=None):
    ''' Execute a task and save its results

    Args:
       task (:obj:`Task`): task
       variables (:obj:`list` of :obj:`Variable`): variables that should be recorded
       log (:obj:`TaskLog`, optional): log for the task

    Returns:
        :obj:`tuple`:

            :obj:`VariableResults`: results of variables
            :obj:`TaskLog`: log
    '''
    config = get_config()

    log = log or TaskLog()

    model = task.model
    sim = task.simulation

    if config.VALIDATE_SEDML:
        raise_errors_warnings(validation.validate_task(task),
                              error_summary='Task `{}` is invalid.'.format(task.id))
        raise_errors_warnings(validation.validate_model_language(model.language, ModelLanguage.SBML),
                              error_summary='Language for model `{}` is not supported.'.format(model.id))
        raise_errors_warnings(validation.validate_model_change_types(model.changes, ()),
                              error_summary='Changes for model `{}` are not supported.'.format(model.id))
        raise_errors_warnings(*validation.validate_model_changes(task.model),
                              error_summary='Changes for model `{}` are invalid.'.format(model.id))
        raise_errors_warnings(validation.validate_simulation_type(sim, (UniformTimeCourseSimulation, )),
                              error_summary='{} `{}` is not supported.'.format(sim.__class__.__name__, sim.id))
        raise_errors_warnings(*validation.validate_simulation(sim),
                              error_summary='Simulation `{}` is invalid.'.format(sim.id))
        raise_errors_warnings(*validation.validate_data_generator_variables(variables),
                              error_summary='Data generator variables for task `{}` are invalid.'.format(task.id))

    target_x_paths_to_sbml_ids = validation.validate_variable_xpaths(variables, model.source, attr='id')

    if config.VALIDATE_SEDML_MODELS:
        raise_errors_warnings(*validation.validate_model(model, [], working_dir='.'),
                              error_summary='Model `{}` is invalid.'.format(model.id),
                              warning_summary='Model `{}` may be invalid.'.format(model.id))

    # validate time course
    if sim.initial_time != 0:
        msg = 'Initial time must be zero, not `{}`.'.format(sim.initial_time)
        raise NotImplementedError(msg)

    # determine the simulation algorithm
    algorithm_substitution_policy = get_algorithm_substitution_policy()
    exec_kisao_id = get_preferred_substitute_algorithm_by_ids(
        sim.algorithm.kisao_id, KISAO_ALGORITHMS_MAP.keys(),
        substitution_policy=algorithm_substitution_policy)
    if exec_kisao_id == sim.algorithm.kisao_id:
        if (
            ALGORITHM_SUBSTITUTION_POLICY_LEVELS[algorithm_substitution_policy]
            > ALGORITHM_SUBSTITUTION_POLICY_LEVELS[AlgorithmSubstitutionPolicy.NONE]
        ):
            changes = []
            unsupported_changes = []
            for change in sim.algorithm.changes:
                if change.kisao_id in ['KISAO_0000594', 'KISAO_0000483']:
                    changes.append(change)
                else:
                    unsupported_changes.append(change.kisao_id)
            if unsupported_changes:
                warn('{} unsuported algorithm parameters were ignored:\n  {}'.format(
                    len(unsupported_changes), '\n  '.join(sorted(unsupported_changes))),
                    BioSimulatorsWarning)
        else:
            changes = sim.algorithm.changes
    else:
        changes = []
    algorithm = Algorithm(kisao_id=exec_kisao_id, changes=changes)

    # determine the simulation method and its parameters
    integrator, time_step = get_integrator(algorithm)

    if time_step is None:
        time_step = (sim.output_end_time - sim.output_start_time) / sim.number_of_steps / 1e2
    use_lazy_newton_method = 0

    number_of_steps = (sim.output_end_time - sim.initial_time) / time_step
    if abs(number_of_steps - round(number_of_steps)) > 1e-8:
        msg = (
            'Time course and time step must specify a positive integer number of steps, not `{}`.'
            '\n  Initial time: {}'
            '\n  Output start time: {}'
            '\n  Output end time: {}'
            '\n  Number of steps: {}'
            '\n  Time step: {}'
        ).format(
            number_of_steps,
            sim.initial_time,
            sim.output_start_time,
            sim.output_end_time,
            sim.number_of_steps,
            time_step,
        )
        raise ValueError(msg)

    print_interval = (sim.output_end_time - sim.output_start_time) / sim.number_of_steps / time_step
    if print_interval <= 0 or abs(print_interval - round(print_interval)) > 1e-8:
        msg = (
            'Time course and time step must specify a positive integer number of intervals, not `{}`.'
            '\n  Initial time: {}'
            '\n  Output start time: {}'
            '\n  Output end time: {}'
            '\n  Number of steps: {}'
            '\n  Time step: {}'
        ).format(
            print_interval,
            sim.initial_time,
            sim.output_start_time,
            sim.output_end_time,
            sim.number_of_steps,
            time_step,
        )
        raise ValueError(msg)
    print_interval = round(print_interval)

    print_amount = 0

    # execute the simulation
    results = libsbmlsim.simulateSBMLFromFile(model.source,
                                              sim.output_end_time, time_step,
                                              print_interval, print_amount,
                                              integrator, use_lazy_newton_method)

    if results.isError():
        raise ValueError(results.error_message)

    # read results through CSV because the ``get*ValueAtIndex`` functions have bugs
    fid, filename = tempfile.mkstemp(suffix='.csv')
    os.close(fid)
    libsbmlsim.write_csv(results, filename)
    results_df = pandas.read_csv(filename)
    os.remove(filename)

    # extract results
    unsupported_symbols = []
    unsupported_targets = []
    variable_results = VariableResults()
    for variable in variables:
        if variable.symbol:
            if variable.symbol == Symbol.time.value:
                variable_result = results_df['time']
            else:
                variable_result = None
                unsupported_symbols.append((variable.id, variable.symbol))

        else:
            sbml_id = target_x_paths_to_sbml_ids[variable.target]
            if sbml_id in results_df:
                variable_result = results_df[sbml_id]
            else:
                variable_result = None
                unsupported_targets.append((variable.id, variable.target))

        if variable_result is not None:
            if exec_kisao_id in ['KISAO_0000086', 'KISAO_0000321']:
                variable_results[variable.id] = variable_result[-(sim.number_of_steps*print_interval + 1)::print_interval].to_numpy()

            else:
                variable_results[variable.id] = variable_result[-(sim.number_of_steps + 1):].to_numpy()

    if unsupported_symbols:
        msg = '{} variables involve unsupported symbols:\n  {}\n\nThe following symbols are supported:\n  {}'.format(
            len(unsupported_symbols),
            '\n  '.join('{}: {}'.format(id, symbol) for id, symbol in sorted(unsupported_symbols)),
            '\n  '.join([Symbol.time.value]))
        raise NotImplementedError(msg)

    if unsupported_targets:
        supported_targets = []

        for i_compartment in range(results.getNumOfCompartments()):
            id = results.getCompartmentNameAtIndex(i_compartment)
            supported_targets.append("/sbml:sbml/sbml:model/sbml:listOfCompartments/sbml:compartment[@id='{}']".format(id))

        for i_parameter in range(results.getNumOfParameters()):
            id = results.getParameterNameAtIndex(i_parameter)
            supported_targets.append("/sbml:sbml/sbml:model/sbml:listOfParameters/sbml:parameter[@id='{}']".format(id))

        for i_species in range(results.getNumOfSpecies()):
            id = results.getSpeciesNameAtIndex(i_species)
            supported_targets.append("/sbml:sbml/sbml:model/sbml:listOfSpecies/sbml:species[@id='{}']".format(id))

        msg = '{} variables involve unsupported targets:\n  {}\n\nThe following targets are supported:\n  {}'.format(
            len(unsupported_targets),
            '\n  '.join('{}: {}'.format(id, target) for id, target in sorted(unsupported_targets)),
            '\n  '.join(supported_targets))
        raise NotImplementedError(msg)

    # log action
    log.algorithm = exec_kisao_id
    log.simulator_details = {
        'method': "simulateSBMLFromFile",
        'arguments': {
            'sim_time': sim.output_end_time,
            'dt': time_step,
            'print_interval': print_interval,
            'print_amount': print_amount,
            'method': integrator,
            'use_lazy_method': use_lazy_newton_method,
        },
    }

    # return results and log
    return variable_results, log
