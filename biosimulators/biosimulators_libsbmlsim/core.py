""" Methods for executing SED tasks in COMBINE archives and saving their outputs

:Author: Jonathan Karr <karr@mssm.edu>
:Date: 2021-03-27
:Copyright: 2021, Center for Reproducible Biomedical Modeling
:License: MIT
"""

from .data_model import get_integrator
from biosimulators_utils.combine.exec import exec_sedml_docs_in_archive
from biosimulators_utils.log.data_model import CombineArchiveLog, TaskLog  # noqa: F401
from biosimulators_utils.plot.data_model import PlotFormat  # noqa: F401
from biosimulators_utils.report.data_model import ReportFormat, VariableResults  # noqa: F401
from biosimulators_utils.sedml.data_model import (Task, ModelLanguage, UniformTimeCourseSimulation,  # noqa: F401
                                                  Variable, Symbol)
from biosimulators_utils.sedml import validation
from biosimulators_utils.sedml.exec import exec_sed_doc
from biosimulators_utils.utils.core import raise_errors_warnings
import functools
import libsbmlsim


__all__ = ['exec_sedml_docs_in_combine_archive', 'exec_sed_task']


def exec_sedml_docs_in_combine_archive(archive_filename, out_dir,
                                       report_formats=None, plot_formats=None,
                                       bundle_outputs=None, keep_individual_outputs=None):
    """ Execute the SED tasks defined in a COMBINE/OMEX archive and save the outputs

    Args:
        archive_filename (:obj:`str`): path to COMBINE/OMEX archive
        out_dir (:obj:`str`): path to store the outputs of the archive

            * CSV: directory in which to save outputs to files
              ``{ out_dir }/{ relative-path-to-SED-ML-file-within-archive }/{ report.id }.csv``
            * HDF5: directory in which to save a single HDF5 file (``{ out_dir }/reports.h5``),
              with reports at keys ``{ relative-path-to-SED-ML-file-within-archive }/{ report.id }`` within the HDF5 file

        report_formats (:obj:`list` of :obj:`ReportFormat`, optional): report format (e.g., csv or h5)
        plot_formats (:obj:`list` of :obj:`PlotFormat`, optional): report format (e.g., pdf)
        bundle_outputs (:obj:`bool`, optional): if :obj:`True`, bundle outputs into archives for reports and plots
        keep_individual_outputs (:obj:`bool`, optional): if :obj:`True`, keep individual output files

    Returns:
        :obj:`CombineArchiveLog`: log
    """
    sed_doc_executer = functools.partial(exec_sed_doc, exec_sed_task)
    return exec_sedml_docs_in_archive(sed_doc_executer, archive_filename, out_dir,
                                      apply_xml_model_changes=True,
                                      report_formats=report_formats,
                                      plot_formats=plot_formats,
                                      bundle_outputs=bundle_outputs,
                                      keep_individual_outputs=keep_individual_outputs)


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
    log = log or TaskLog()

    raise_errors_warnings(validation.validate_model_language(task.model.language, ModelLanguage.SBML),
                          error_summary='Language for model `{}` is not supported.'.format(model.id))
    raise_errors_warnings(validation.validate_model_change_types(task.model.changes, ()),
                          error_summary='Changes for model `{}` are not supported.'.format(model.id))
    raise_errors_warnings(validation.validate_simulation_type(task.simulation, (UniformTimeCourseSimulation, )),
                          error_summary='{} `{}` is not supported.'.format(sim.__class__.__name__, sim.id))
    target_x_paths_to_sbml_ids = validation.validate_variable_xpaths(variables, task.model.source, attr='id')

    # validate time course
    if task.simulation.initial_time != 0:
        raise NotImplementedError('Initial time must be zero')

    # execute the simulation
    integrator = get_integrator(task.simulation.algorithm)
    args = [task.model,
            task.simulation.output_end_time, task.simulation.output_end_time / task.simulation.number_of_points,
            1, 0, integrator, 0]
    results = libsbmlsim.simulateSBMLFromFile(*args)

    # extract results
    variable_results = VariableResults()
    for variable in variables:
        # TODO: extract results from `results` and convert to numpy array
        # TODO: throw errors in the simulation request an invalid/unsupported output

        if variable.symbol:
            if variable.symbol == Symbol.time:
                i_result = todo_index_of_time_in_results  # index in results for time
            else:
                raise NotImplementedError('Symbol {} is not supported'.format(variable.symbol))

        else:
            sbml_id = target_x_paths_to_sbml_ids[variable.target]
            i_result = todo_index_of_sbml_id_in_results  # index in results for SBML id
            if not i_result:
                raise NotImplementedError('Variable with target {} is not supported'.format(variable.target))

        variable_results[variable.id] = results[i_result][-(task.simulation.number_of_points + 1):]

    # log action
    log.algorithm = task.simulation.algorithm.kisao_id
    log.simulator_details = {
        'method': "simulateSBMLFromFile",
        'arguments': args,
    }

    # return results and log
    return variable_results, log
