from biosimulators_libsbmlsim import data_model
from biosimulators_utils.sedml.data_model import Algorithm, AlgorithmParameterChange
import json
import libsbmlsim
import os
import unittest


class DataModelTestCase(unittest.TestCase):
    def test_get_integrator_parameters(self):
        alg = Algorithm(
            kisao_id='KISAO_0000279',
            changes=[
                AlgorithmParameterChange(kisao_id='KISAO_0000594', new_value='2'),
                AlgorithmParameterChange(kisao_id='KISAO_0000483', new_value='0.01'),
            ],
        )
        self.assertEqual(data_model.get_integrator_parameters(alg, [1, 2, 3, 4]),
                         (2, 0.01))

        alg.changes[0].new_value = '4.1'
        with self.assertRaisesRegex(ValueError, 'must be an integer between'):
            data_model.get_integrator_parameters(alg, [1, 2, 3, 4])

        alg.changes[0].new_value = '5'
        with self.assertRaisesRegex(NotImplementedError, 'does not support order'):
            data_model.get_integrator_parameters(alg, [1, 2, 3, 4])

        alg.changes = alg.changes[1:]
        self.assertEqual(data_model.get_integrator_parameters(alg, [1, 2, 3, 4]),
                         (None, 0.01))

        alg.changes = []
        self.assertEqual(data_model.get_integrator_parameters(alg, [1, 2, 3, 4]),
                         (None, None))

        alg = Algorithm(
            kisao_id='KISAO_0000030',
            changes=[
                AlgorithmParameterChange(kisao_id='KISAO_0000483', new_value='0.01'),
            ],
        )
        self.assertEqual(data_model.get_integrator_parameters(alg, None),
                         (None, 0.01))

        alg.changes[0].new_value = 'abc'
        with self.assertRaisesRegex(ValueError, 'must be a positive'):
            data_model.get_integrator_parameters(alg, None)

        alg.changes[0].new_value = '-5'
        with self.assertRaisesRegex(ValueError, 'must be a positive'):
            data_model.get_integrator_parameters(alg, None)

        alg.changes[0].kisao_id = 'KISAO_0000594'
        alg.changes[0].new_value = '1'
        with self.assertRaisesRegex(NotImplementedError, 'does not support parameter'):
            data_model.get_integrator_parameters(alg, None)

    def test_get_integrator(self):
        alg = Algorithm(
            kisao_id='KISAO_0000279',
            changes=[
                AlgorithmParameterChange(kisao_id='KISAO_0000594', new_value='2'),
                AlgorithmParameterChange(kisao_id='KISAO_0000483', new_value='0.01'),
            ],
        )
        self.assertEqual(data_model.get_integrator(alg), (libsbmlsim.MTHD_ADAMS_BASHFORTH_2, 0.01))

        alg.changes = alg.changes[1:]
        self.assertEqual(data_model.get_integrator(alg), (libsbmlsim.MTHD_ADAMS_BASHFORTH_4, 0.01))

        alg.changes = []
        self.assertEqual(data_model.get_integrator(alg), (libsbmlsim.MTHD_ADAMS_BASHFORTH_4, None))

        alg.kisao_id = 'KISAO_0000030'
        self.assertEqual(data_model.get_integrator(alg), (libsbmlsim.MTHD_EULER, None))

        for alg_props in data_model.KISAO_ALGORITHMS_MAP.values():
            alg = Algorithm(kisao_id=alg_props['id'])
            alg.changes.append(AlgorithmParameterChange(kisao_id='KISAO_0000483', new_value='0.01'))
            method, time_step = data_model.get_integrator(alg)
            self.assertIsInstance(method, int)
            self.assertEqual(time_step, 0.01)

        alg.kisao_id = 'KISAO_0000086'
        alg.changes.append(AlgorithmParameterChange(kisao_id='KISAO_0000483', new_value='0.01'))
        data_model.get_integrator(alg)

        alg.kisao_id = 'KISAO_0000019'
        with self.assertRaisesRegex(NotImplementedError, 'supports the following algorithms'):
            data_model.get_integrator(alg)

    def test_consistent_with_specs(self):
        with open(os.path.join(os.path.dirname(__file__), '..', '..', 'biosimulators.json'), 'r') as file:
            specs = json.load(file)

        self.assertEqual(set(data_model.KISAO_ALGORITHMS_MAP.keys()),
                         set(alg_specs['kisaoId']['id'] for alg_specs in specs['algorithms']))

        for alg_specs in specs['algorithms']:
            alg_props = data_model.KISAO_ALGORITHMS_MAP[alg_specs['kisaoId']['id']]

            param_kisao_ids = set()
            param_kisao_ids.add('KISAO_0000483')
            if len(alg_props['orders']) > 1:
                param_kisao_ids.add('KISAO_0000594')
            self.assertEqual(param_kisao_ids,
                             set(param_specs['kisaoId']['id'] for param_specs in alg_specs['parameters']))

            for param_specs in alg_specs['parameters']:
                if param_specs['kisaoId']['id'] == 'KISAO_0000594':
                    self.assertEqual(param_specs['type'], 'integer')

                    range = list(alg_props['orders'].keys())
                    range.remove(None)
                    range.sort()

                    self.assertEqual(param_specs['recommendedRange'], [str(range[0]), str(range[-1])])
                    self.assertEqual(param_specs['value'], str(range[-1]))

                elif param_specs['kisaoId']['id'] == 'KISAO_0000483':
                    self.assertEqual(param_specs['type'], 'float')
                    self.assertEqual(param_specs['recommendedRange'], None)
                    self.assertEqual(param_specs['value'], None)
