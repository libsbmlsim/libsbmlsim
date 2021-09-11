""" Methods for executing SED tasks in COMBINE archives and saving their outputs

:Author: Jonathan Karr <karr@mssm.edu>
:Date: 2021-03-27
:Copyright: 2021, Center for Reproducible Biomedical Modeling
:License: MIT
"""

from .data_model import KISAO_ALGORITHMS_MAP, get_integrator
from biosimulators_utils.combine.exec import exec_sedml_docs_in_archive
from biosimulators_utils.config import get_config, Config  # noqa: F401
from biosimulators_utils.log.data_model import CombineArchiveLog, TaskLog, StandardOutputErrorCapturerLevel  # noqa: F401
from biosimulators_utils.viz.data_model import VizFormat  # noqa: F401
from biosimulators_utils.report.data_model import ReportFormat, VariableResults, SedDocumentResults  # noqa: F401
from biosimulators_utils.sedml.data_model import (Task, ModelLanguage, ModelAttributeChange, UniformTimeCourseSimulation,  # noqa: F401
                                                  Algorithm, Variable, Symbol)
from biosimulators_utils.sedml import validation
from biosimulators_utils.sedml.exec import exec_sed_doc as base_exec_sed_doc
from biosimulators_utils.sedml.utils import apply_changes_to_xml_model
from biosimulators_utils.simulator.utils import get_algorithm_substitution_policy
from biosimulators_utils.utils.core import raise_errors_warnings
from biosimulators_utils.warnings import warn, BioSimulatorsWarning
from kisao.data_model import AlgorithmSubstitutionPolicy, ALGORITHM_SUBSTITUTION_POLICY_LEVELS
from kisao.utils import get_preferred_substitute_algorithm_by_ids
import copy
import libsbmlsim
import lxml.etree
import os
import pandas
import tempfile


__all__ = ['exec_sedml_docs_in_combine_archive', 'exec_sed_doc', 'exec_sed_task', 'preprocess_sed_task']


def exec_sedml_docs_in_combine_archive(archive_filename, out_dir, config=None):
    """ Execute the SED tasks defined in a COMBINE/OMEX archive and save the outputs

    Args:
        archive_filename (:obj:`str`): path to COMBINE/OMEX archive
        out_dir (:obj:`str`): path to store the outputs of the archive

            * CSV: directory in which to save outputs to files
              ``{ out_dir }/{ relative-path-to-SED-ML-file-within-archive }/{ report.id }.csv``
            * HDF5: directory in which to save a single HDF5 file (``{ out_dir }/reports.h5``),
              with reports at keys ``{ relative-path-to-SED-ML-file-within-archive }/{ report.id }`` within the HDF5 file

        config (:obj:`Config`, optional): BioSimulators common configuration

    Returns:
        :obj:`tuple`:

            * :obj:`SedDocumentResults`: results
            * :obj:`CombineArchiveLog`: log
    """
    return exec_sedml_docs_in_archive(exec_sed_doc, archive_filename, out_dir,
                                      apply_xml_model_changes=True,
                                      config=config)


def exec_sed_doc(doc, working_dir, base_out_path, rel_out_path=None,
                 apply_xml_model_changes=True,
                 log=None, indent=0, pretty_print_modified_xml_models=False,
                 log_level=StandardOutputErrorCapturerLevel.c, config=None):
    """ Execute the tasks specified in a SED document and generate the specified outputs

    Args:
        doc (:obj:`SedDocument` or :obj:`str`): SED document or a path to SED-ML file which defines a SED document
        working_dir (:obj:`str`): working directory of the SED document (path relative to which models are located)

        base_out_path (:obj:`str`): path to store the outputs

            * CSV: directory in which to save outputs to files
              ``{base_out_path}/{rel_out_path}/{report.id}.csv``
            * HDF5: directory in which to save a single HDF5 file (``{base_out_path}/reports.h5``),
              with reports at keys ``{rel_out_path}/{report.id}`` within the HDF5 file

        rel_out_path (:obj:`str`, optional): path relative to :obj:`base_out_path` to store the outputs
        apply_xml_model_changes (:obj:`bool`, optional): if :obj:`True`, apply any model changes specified in the SED-ML file before
            calling :obj:`task_executer`.
        log (:obj:`SedDocumentLog`, optional): log of the document
        indent (:obj:`int`, optional): degree to indent status messages
        pretty_print_modified_xml_models (:obj:`bool`, optional): if :obj:`True`, pretty print modified XML models
        log_level (:obj:`StandardOutputErrorCapturerLevel`, optional): level at which to log output
        config (:obj:`Config`, optional): BioSimulators common configuration
        simulator_config (:obj:`SimulatorConfig`, optional): tellurium configuration

    Returns:
        :obj:`tuple`:

            * :obj:`ReportResults`: results of each report
            * :obj:`SedDocumentLog`: log of the document
    """
    return base_exec_sed_doc(exec_sed_task, doc, working_dir, base_out_path,
                             rel_out_path=rel_out_path,
                             apply_xml_model_changes=apply_xml_model_changes,
                             log=log,
                             indent=indent,
                             pretty_print_modified_xml_models=pretty_print_modified_xml_models,
                             log_level=log_level,
                             config=config)


def exec_sed_task(task, variables, preprocessed_task=None, log=None, config=None):
    ''' Execute a task and save its results

    Args:
        task (:obj:`Task`): task
        variables (:obj:`list` of :obj:`Variable`): variables that should be recorded
        preprocessed_task (:obj:`object`, optional): preprocessed information about the task, including possible
            model changes and variables. This can be used to avoid repeatedly executing the same initialization
            for repeated calls to this method.
        log (:obj:`TaskLog`, optional): log for the task
        config (:obj:`Config`, optional): BioSimulators common configuration

    Returns:
        :obj:`tuple`:

            :obj:`VariableResults`: results of variables
            :obj:`TaskLog`: log
    '''
    config = config or get_config()

    if config.LOG and not log:
        log = TaskLog()

    if preprocessed_task is None:
        preprocessed_task = preprocess_sed_task(task, variables, config=config)

    model = task.model
    sim = task.simulation

    # change model
    if model.changes:
        raise_errors_warnings(validation.validate_model_change_types(model.changes, (ModelAttributeChange,)),
                              error_summary='Changes for model `{}` are not supported.'.format(model.id))

        model_etree = preprocessed_task['model']['etree']

        model = copy.deepcopy(model)
        for change in model.changes:
            change.new_value = str(change.new_value)

        apply_changes_to_xml_model(model, model_etree, sed_doc=None, working_dir=None)

        model_file, model_filename = tempfile.mkstemp(suffix='.xml')
        os.close(model_file)

        model_etree.write(model_filename,
                          xml_declaration=True,
                          encoding="utf-8",
                          standalone=False,
                          pretty_print=False)
    else:
        model_filename = model.source

    # validate time course
    if sim.initial_time != 0:
        msg = 'Initial time must be zero, not `{}`.'.format(sim.initial_time)

        if model.changes:
            os.remove(model_filename)

        raise NotImplementedError(msg)

    # determine the number of simulation steps
    time_step = preprocessed_task['simulation']['time_step']

    if time_step is None:
        time_step = (sim.output_end_time - sim.output_start_time) / sim.number_of_steps / 1e2

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

        if model.changes:
            os.remove(model_filename)

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

        if model.changes:
            os.remove(model_filename)

        raise ValueError(msg)
    print_interval = round(print_interval)

    # execute the simulation
    results = libsbmlsim.simulateSBMLFromFile(model_filename,
                                              sim.output_end_time,
                                              time_step,
                                              print_interval,
                                              preprocessed_task['simulation']['print_amount'],
                                              preprocessed_task['simulation']['integrator'],
                                              preprocessed_task['simulation']['use_lazy_newton_method'])

    if model.changes:
        os.remove(model_filename)

    if results.isError():
        raise ValueError(results.error_message)

    # read results through CSV because the ``get*ValueAtIndex`` functions have bugs
    fid, filename = tempfile.mkstemp(suffix='.csv')
    os.close(fid)
    libsbmlsim.write_csv(results, filename)
    results_df = pandas.read_csv(filename)
    os.remove(filename)

    # extract results
    xpath_sbml_id_map = preprocessed_task['model']['xpath_sbml_id_map']
    variable_results = VariableResults()
    unsupported_symbols = []
    unsupported_targets = []
    for variable in variables:
        if variable.symbol:
            if variable.symbol == Symbol.time.value:
                variable_result = results_df['time']
            else:
                variable_result = None
                unsupported_symbols.append((variable.id, variable.symbol))

        else:
            sbml_id = xpath_sbml_id_map[variable.target]
            if sbml_id in results_df:
                variable_result = results_df[sbml_id]
            else:
                variable_result = None
                unsupported_targets.append((variable.id, variable.target))

        if variable_result is not None:
            if preprocessed_task['simulation']['algorithm_kisao_id'] in ['KISAO_0000086', 'KISAO_0000321']:
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
    if config.LOG:
        log.algorithm = preprocessed_task['simulation']['algorithm_kisao_id'],
        log.simulator_details = {
            'method': "simulateSBMLFromFile",
            'arguments': {
                'sim_time': sim.output_end_time,
                'dt': time_step,
                'print_interval': print_interval,
                'print_amount': preprocessed_task['simulation']['print_amount'],
                'method': preprocessed_task['simulation']['integrator'],
                'use_lazy_method': preprocessed_task['simulation']['use_lazy_newton_method'],
            },
        }

    # return results and log
    return variable_results, log


def preprocess_sed_task(task, variables, config=None):
    """ Preprocess a SED task, including its possible model changes and variables. This is useful for avoiding
    repeatedly initializing tasks on repeated calls of :obj:`exec_sed_task`.

    Args:
        task (:obj:`Task`): task
        variables (:obj:`list` of :obj:`Variable`): variables that should be recorded
        config (:obj:`Config`, optional): BioSimulators common configuration

    Returns:
        :obj:`object`: preprocessed information about the task
    """
    config = config or get_config()

    model = task.model
    sim = task.simulation

    # validate model and simulation
    if config.VALIDATE_SEDML:
        raise_errors_warnings(validation.validate_task(task),
                              error_summary='Task `{}` is invalid.'.format(task.id))
        raise_errors_warnings(validation.validate_model_language(model.language, ModelLanguage.SBML),
                              error_summary='Language for model `{}` is not supported.'.format(model.id))
        raise_errors_warnings(validation.validate_model_change_types(model.changes, (ModelAttributeChange,)),
                              error_summary='Changes for model `{}` are not supported.'.format(model.id))
        raise_errors_warnings(*validation.validate_model_changes(task.model),
                              error_summary='Changes for model `{}` are invalid.'.format(model.id))
        raise_errors_warnings(validation.validate_simulation_type(sim, (UniformTimeCourseSimulation, )),
                              error_summary='{} `{}` is not supported.'.format(sim.__class__.__name__, sim.id))
        raise_errors_warnings(*validation.validate_simulation(sim),
                              error_summary='Simulation `{}` is invalid.'.format(sim.id))
        raise_errors_warnings(*validation.validate_data_generator_variables(variables),
                              error_summary='Data generator variables for task `{}` are invalid.'.format(task.id))

    if not os.path.isfile(model.source):
        raise FileNotFoundError('Model source `{}` is not a file.'.format(model.source))
    model_etree = lxml.etree.parse(model.source)
    xpath_sbml_id_map = validation.validate_target_xpaths(variables, model_etree, attr='id')

    if config.VALIDATE_SEDML_MODELS:
        raise_errors_warnings(*validation.validate_model(model, [], working_dir='.'),
                              error_summary='Model `{}` is invalid.'.format(model.id),
                              warning_summary='Model `{}` may be invalid.'.format(model.id))

    # determine the simulation algorithm
    algorithm_substitution_policy = get_algorithm_substitution_policy(config=config)
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
    use_lazy_newton_method = 0
    print_amount = 0

    # return preprocessed task
    return {
        'model': {
            'etree': model_etree,
            'xpath_sbml_id_map': xpath_sbml_id_map,
        },
        'simulation': {
            'algorithm_kisao_id': exec_kisao_id,
            'integrator': integrator,
            'use_lazy_newton_method': use_lazy_newton_method,
            'time_step': time_step,
            'print_amount': print_amount,
        },
    }
