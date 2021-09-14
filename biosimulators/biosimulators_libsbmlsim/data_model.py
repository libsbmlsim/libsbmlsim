""" Data model for mapping KiSAO terms to LibSBMLsim integrators and their settings

:Author: Jonathan Karr <karr@mssm.edu>
:Date: 2021-03-27
:Copyright: 2021, Center for Reproducible Biomedical Modeling
:License: MIT
"""

from biosimulators_utils.sedml.data_model import Algorithm, AlgorithmParameterChange  # noqa: F401
import collections
import libsbmlsim

__all__ = [
    'KISAO_ALGORITHMS_MAP',
    'get_integrator',
    'get_integrator_parameters',
]

KISAO_ALGORITHMS_MAP = collections.OrderedDict([
    ('KISAO_0000086', {
        'id': 'KISAO_0000086',
        'name': 'Fehlberg method',
        'orders': {
            None: libsbmlsim.MTHD_RUNGE_KUTTA_FEHLBERG_5,
        },
        'uses_print_interval': False,
    }),
    ('KISAO_0000030', {
        'id': 'KISAO_0000030',
        'name': 'Forward Euler method',
        'orders': {
            None: libsbmlsim.MTHD_EULER,
        },
        'uses_print_interval': True,
    }),
    ('KISAO_0000279', {
        'id': 'KISAO_0000279',
        'name': 'Adams-Bashforth method',
        'orders': {
            1: libsbmlsim.MTHD_ADAMS_BASHFORTH_1,
            2: libsbmlsim.MTHD_ADAMS_BASHFORTH_2,
            3: libsbmlsim.MTHD_ADAMS_BASHFORTH_3,
            4: libsbmlsim.MTHD_ADAMS_BASHFORTH_4,
            None: libsbmlsim.MTHD_ADAMS_BASHFORTH_4,
        },
        'uses_print_interval': True,
    }),
    ('KISAO_0000032', {
        'id': 'KISAO_0000032',
        'name': 'Explicit fourth-order Runge-Kutta method',
        'orders': {
            None: libsbmlsim.MTHD_RUNGE_KUTTA,
        },
        'uses_print_interval': True,
    }),
    ('KISAO_0000321', {
        'id': 'KISAO_0000321',
        'name': 'Cash-Karp method',
        'orders': {
            None: libsbmlsim.MTHD_CASH_KARP,
        },
        'uses_print_interval': False,
    }),
    ('KISAO_0000031', {
        'id': 'KISAO_0000031',
        'name': 'Backward Euler method',
        'orders': {
            None: libsbmlsim.MTHD_BACKWARD_EULER,
        },
        'uses_print_interval': True,
    }),
    ('KISAO_0000309', {
        'id': 'KISAO_0000309',
        'name': 'Crank-Nicolson method',
        'orders': {
            None: libsbmlsim.MTHD_CRANK_NICOLSON,
        },
        'uses_print_interval': True,
    }),
    ('KISAO_0000280', {
        'id': 'KISAO_0000280',
        'name': 'Adams Moulton method',
        'orders': {
            2: libsbmlsim.MTHD_ADAMS_MOULTON_2,
            3: libsbmlsim.MTHD_ADAMS_MOULTON_3,
            4: libsbmlsim.MTHD_ADAMS_MOULTON_4,
            None: libsbmlsim.MTHD_ADAMS_MOULTON_4,
        },
        'uses_print_interval': True,
    }),
    ('KISAO_0000288', {
        'id': 'KISAO_0000288',
        'name': 'Backward difference method',
        'orders': {
            2: libsbmlsim.MTHD_BACKWARD_DIFFERENCE_2,
            3: libsbmlsim.MTHD_BACKWARD_DIFFERENCE_3,
            4: libsbmlsim.MTHD_BACKWARD_DIFFERENCE_4,
            None: libsbmlsim.MTHD_BACKWARD_DIFFERENCE_4,
        },
        'uses_print_interval': True,
    }),
])


def get_integrator(algorithm):
    """ Get the LiSBMLsim integrator and its parameter for a SED-ML/KiSAO algorithm

    Args:
        algorithm (:obj:`Algorithm`): SED-ML algorithm

    Returns:
        :obj:`tuple`:

            :obj:`int`: id of LibSBMLsim integrator
            :obj:`float`: time step
    """
    method_group = KISAO_ALGORITHMS_MAP.get(algorithm.kisao_id, None)
    if not method_group:
        msg = 'Algorithm {} is not supported. libSBMLSim supports the following algorithms:\n  {}'.format(
            algorithm.kisao_id,
            '\n  '.join(sorted('{}: {}'.format(alg['name'], alg['id']) for alg in KISAO_ALGORITHMS_MAP.values()))
        )
        raise NotImplementedError(msg)
    supported_orders = list(method_group['orders'].keys())
    supported_orders.remove(None)
    order, dt = get_integrator_parameters(algorithm,
                                          supported_orders=supported_orders)
    method = method_group['orders'][order]
    return (method, dt)


def get_integrator_parameters(algorithm, supported_orders=None):
    """ Get the order and time step of a LiSBMLsim integrator for a SED-ML/KiSAO algorithm

    Args:
        algorithm (:obj:`Algorithm`): SED-ML algorithm
        supported_orders (:obj:`list` of :obj:`int`, optional): orders supported by the SED-ML algorithm

    Returns:
        :obj:`tuple`:

            * :obj:`int`: order
            * :obj:`float`: time step
    """
    order = None
    time_step = None
    for change in algorithm.changes:
        if supported_orders and change.kisao_id == 'KISAO_0000594':
            try:
                order = int(change.new_value)
            except ValueError:
                raise ValueError('The order of algorithm {} must be an integer between {} and {}.'.format(
                    algorithm.kisao_id, supported_orders[0], supported_orders[-1]))

            if order not in supported_orders:
                raise NotImplementedError('Algorithm {} does not support order {} (KISAO_0000594).',
                                          algorithm.kisao_id, order)

        elif change.kisao_id == 'KISAO_0000483':
            try:
                time_step = float(change.new_value)
            except ValueError:
                raise ValueError('The time step must be a positive float, not `{}`.'.format(
                    change.new_value))

            if time_step <= 0:
                raise ValueError('The time step must be a positive float, not `{}`.'.format(
                    change.new_value))

        else:
            raise NotImplementedError('Algorithm {} does not support parameter {}.'.format(
                algorithm.kisao_id, change.kisao_id))

    return (order, time_step)
