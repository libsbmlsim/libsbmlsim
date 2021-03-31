""" Data model for mapping KiSAO terms to LibSBMLsim integrators and their settings

:Author: Jonathan Karr <karr@mssm.edu>
:Date: 2021-03-27
:Copyright: 2021, Center for Reproducible Biomedical Modeling
:License: MIT
"""

from biosimulators_utils.sedml.data_model import Algorithm, AlgorithmParameterChange  # noqa: F401
import libsbmlsim

__all__ = [
    'get_integrator',
    'get_order',
]


def get_integrator(algorithm):
    """ Get the LiSBMLsim integrator for a SED-ML/KiSAO algorithm

    Args:
        algorithm (:obj:`Algorithm`): SED-ML algorithm

    Returns:
        :obj:`int`: id of LibSBMLsim integrator
    """
    if algorithm.kisao_id == 'KISAO_0000030':
        get_order(algorithm)
        return libsbmlsim.MTHD_EULER

    elif algorithm.kisao_id == 'KISAO_0000279':
        order = get_order(algorithm, supported_orders=[1, 2, 3, 4])
        if order == 1:
            return libsbmlsim.MTHD_ADAMS_BASHFORTH_1
        elif order == 2:
            return libsbmlsim.MTHD_ADAMS_BASHFORTH_2
        elif order == 3:
            return libsbmlsim.MTHD_ADAMS_BASHFORTH_3
        elif order == 4:
            return libsbmlsim.MTHD_ADAMS_BASHFORTH_4

    elif algorithm.kisao_id == 'KISAO_0000032':
        get_order(algorithm)
        return libsbmlsim.MTHD_RUNGE_KUTTA

    elif algorithm.kisao_id == 'KISAO_0000086':
        get_order(algorithm)
        return libsbmlsim.MTHD_RUNGE_KUTTA_FEHLBERG_5

    elif algorithm.kisao_id == 'KISAO_0000321':
        get_order(algorithm)
        return libsbmlsim.MTHD_CASH_KARP

    elif algorithm.kisao_id == 'KISAO_0000031':
        get_order(algorithm)
        return libsbmlsim.MTHD_BACKWARD_EULER

    elif algorithm.kisao_id == 'KISAO_0000309':
        get_order(algorithm)
        return libsbmlsim.MTHD_CRANK_NICOLSON

    elif algorithm.kisao_id == 'KISAO_0000280':
        order = get_order(algorithm, supported_orders=[2, 3, 4])
        if order == 2:
            return libsbmlsim.MTHD_ADAMS_MOULTON_2
        elif order == 3:
            return libsbmlsim.MTHD_ADAMS_MOULTON_3
        elif order == 4:
            return libsbmlsim.MTHD_ADAMS_MOULTON_4

    elif algorithm.kisao_id == 'KISAO_0000288':
        order = get_order(algorithm, supported_orders=[2, 3, 4])
        if order == 2:
            return libsbmlsim.MTHD_BACKWARD_DIFFERENCE_2
        elif order == 3:
            return libsbmlsim.MTHD_BACKWARD_DIFFERENCE_3
        elif order == 4:
            return libsbmlsim.MTHD_BACKWARD_DIFFERENCE_4

    else:
        raise NotImplementedError('Algorithm {} is not suppoorted'.format(algorithm.kisao_id))


def get_order(algorithm, supported_orders=None):
    """ Get the desired order of a LiSBMLsim integrator for a SED-ML/KiSAO algorithm

    Args:
        algorithm (:obj:`Algorithm`): SED-ML algorithm
        supported_orders (:obj:`list` of :obj:`int`, optional): orders supported by the SED-ML algorithm

    Returns:
        :obj:`int`: order
    """
    order = None
    for change in algorithm.changes:
        if supported_orders and change.kisao_id == 'KISAO_0000594':
            try:
                order = int(change.new_value)
            except ValueError:
                raise ValueError('The order of algorithm {} must be an integer between {} and {}'.format(
                    algorithm.kisao_id, supported_orders[0], supported_orders[-1]))

            if order not in supported_orders:
                raise NotImplementedError('Algorithm {} does not supported order {}', algorithm.kisao_id, order)

        else:
            raise NotImplementedError('Algorithm {} does not support algorithm {}'.format(
                algorithm.kisao_id, change.kisao_id))

    return order
