import os
import sys
import json
import typing as t

import control as ct
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from numpy import cos, sin, sqrt, tan
from pycomex.experiment import Experiment
from pycomex.util import Skippable
from scipy.optimize import minimize

import labor_regelungstechnik.systems.single_pendulum_nonlinear as single_pendulum_nonlinear
from labor_regelungstechnik.utils import MEASUREMENTS_PATH

# == DATA PARAMETERS ==
MEASUREMENTS_JSON_PATH = os.path.join(MEASUREMENTS_PATH, 'measurements_001.json')

# == SYSTEM PARAMETERS ==
IO_SYSTEM = single_pendulum_nonlinear.io_system

# == EXPERIMENT PARAMETERS ==
BASE_PATH = os.getcwd()
NAMESPACE = 'parameter_optimization'
DEBUG = True
with Skippable(), (e := Experiment(base_path=BASE_PATH, namespace=NAMESPACE, glob=globals())):

    e.info('starting parameter identification...')

    # -- SETTING UP THE SYSTEM --

    # -- LOADING THE MEASUREMENTS --
    with open(MEASUREMENTS_JSON_PATH, mode='r') as file:
        content = file.read()
        measurements: t.List[dict] = json.loads(content)

    # -- OBJECTIVE FUNCTION BASED ON MEASUREMENTS --
    def objective_function(parameters: t.Sequence[float],
                           input_keys: t.Sequence[str] = ('x_const_mess', 'y_const_mess'),
                           input_delays: t.Sequence[float] = (0.3, 0.1),
                           state_keys: t.Sequence[str] = (None, None, 'y_out_mess', None, 'phi_out_mess',
                                                          'x_out_mess'),
                           output_keys: t.Sequence[str] = ('x_out_mess', 'y_out_mess', 'phi_out_mess'),
                           output_weights: t.Sequence[float] = (0.2, 0.2, 1),
                           return_records: bool = False
                           ):
        record_dict = {}

        total_error = 0
        for index, measurement in enumerate(measurements):
            record_dict[index] = {'measured': [], 'simulated': []}

            # We extract the exact same time frame from the model as we do for the measurements
            ts = np.array(measurement['timestamps'])
            record_dict[index]['timestamps'] = ts

            # Then we need to extract the input signals from the measurement
            inputs = []
            for key, delay in zip(input_keys, input_delays):
                values = []
                for t, value in zip(ts, measurement[key]):
                    if t < delay:
                        values.append(0)
                    else:
                        values.append(value)

                inputs.append(np.array(values))

            # Then we also need to extract the starting conditions for the states from the measurements
            initial_conditions = []
            for key in state_keys:
                if key is None:
                    initial_conditions.append(0)
                else:
                    value = measurement[key][3]
                    initial_conditions.append(value)

            initial_conditions[2] = initial_conditions[2]
            initial_conditions[4] = np.radians(initial_conditions[4])

            # We unpack the parameter array into a dict such that the system can make sense of it
            if parameters is None:
                params = {}
            else:
                params = {
                    'm_x': parameters[0],
                    'm_y': parameters[1],
                    'c_varphi': parameters[2],
                }
            # Then we can simulate the model in this time frame
            try:
                _, y = ct.input_output_response(
                    IO_SYSTEM,
                    ts,
                    U=inputs,
                    X0=initial_conditions,
                    solve_ivp_kwargs={
                        #'method': 'LSODA'
                    },
                    params=params
                )
            except RuntimeError:
                return 1000

            # Then we can compare the measurements
            for key, weight, values_simulated in zip(output_keys, output_weights, y):
                values_measured = np.array(measurement[key])
                total_error += weight * np.mean(np.abs(values_measured - values_simulated))

                record_dict[index]['measured'].append(values_measured)
                record_dict[index]['simulated'].append(values_simulated)

        if return_records:
            return total_error, record_dict
        else:
            del record_dict
            return total_error

    # -- PLOTTING THE DEFAULT PARAMETERS --
    error, records_map = objective_function(None, return_records=True)
    e.info(f'mse with default parameters is: {error:.2f}')

    pdf_path = os.path.join(e.path, 'default_parameters.pdf')
    with PdfPages(pdf_path) as pdf:
        for index, records in records_map.items():
            n_rows = 3
            fig, rows = plt.subplots(ncols=1, nrows=n_rows, figsize=(16, 6*n_rows), squeeze=False)
            fig.suptitle(f'measurement {index}\n'
                         f'error: {error:.2f}')

            for row_index, name, limits, row in zip([0, 1, 2],
                                                    ['x', 'l', 'phi'],
                                                    [(-0.1, 3), (-0.1, 1.5), (-15, 15)],
                                                    rows):
                ax = row[0]
                ax.set_title(name)

                ax.set_ylim(limits)
                ax.plot(
                    records['timestamps'],
                    records['measured'][row_index],
                    color='gray',
                    label='measured'
                )
                ax.plot(
                    records['timestamps'],
                    records['simulated'][row_index],
                    color='blue',
                    label='simulated'
                )
                ax.legend()

            pdf.savefig(fig)
            plt.close(fig)

    sys.exit(0)
    # -- PARAMETER OPTIMIZATION --
    e.info('starting parameter optimization...')
    initial_parameters = [35, 3.1, 0.7]
    result = minimize(
        objective_function,
        initial_parameters,
        method='nelder-mead',
        options={
            'maxiter': 10,
            'xatol': 1e-2,
            'disp': True
        }
    )

    # -- PLOTTING OPTIMIZED --
    error, records_map = objective_function(result.x, return_records=True)
    e.info(f'mse with optimized parameters is: {error:.2f}')
    e.info(f'optimized parameters: {result.x}')

    pdf_path = os.path.join(e.path, 'optimized_parameters.pdf')
    with PdfPages(pdf_path) as pdf:
        for index, records in records_map.items():
            n_rows = 3
            fig, rows = plt.subplots(ncols=1, nrows=n_rows, figsize=(16, 6 * n_rows), squeeze=False)
            fig.suptitle(f'measurement {index}\n'
                         f'error: {error:.2f}')

            for row_index, name, row in zip([0, 1, 2], ['x', 'l', 'phi'], rows):
                ax = row[0]
                ax.set_title(name)

                ax.plot(
                    records['timestamps'],
                    records['measured'][row_index],
                    color='gray',
                    label='measured'
                )
                ax.plot(
                    records['timestamps'],
                    records['simulated'][row_index],
                    color='blue',
                    label='simulated'
                )
                ax.legend()

            pdf.savefig(fig)
            plt.close(fig)
