import os
import re
from pprint import pprint, PrettyPrinter
from collections import defaultdict

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from pycomex.experiment import Experiment
from pycomex.util import Skippable
from asammdf import MDF

from labor_regelungstechnik.utils import MEASUREMENTS_PATH

MEASUREMENTS_FILE_NAME = 'messungen_18_11_22.mf4'

BASE_PATH = os.getcwd()
NAMESPACE = 'extract_measurements'
DEBUG = True
with Skippable(), (e := Experiment(base_path=BASE_PATH, namespace=NAMESPACE, glob=globals())):

    measurements_file_path = os.path.join(MEASUREMENTS_PATH, MEASUREMENTS_FILE_NAME)
    e.info(f'starting to extract measurements from recorded file "{measurements_file_path}"')

    # - READING FROM FILE AND PLOTTING

    e.info(f'printing information overview...')
    mdf = MDF(measurements_file_path)
    mdf_info = mdf.info()
    e['info'] = mdf_info
    pretty_printer = PrettyPrinter(indent=4)
    e.info(pretty_printer.pformat(mdf_info))

    signal_names = []
    for value in mdf_info['group 0'].values():
        if isinstance(value, str) and 'name' in value:
            pattern = re.compile(r'name="(.*)"')
            m = re.search(pattern, value)
            if m:
                name = m.group(1)
                signal_names.append(name)

    e.info('plotting all the signals...')
    n_rows = len(signal_names)
    fig, rows = plt.subplots(ncols=1, nrows=n_rows, squeeze=False, figsize=(16, 4 * n_rows))

    signal_map = {}
    for i, signal_name in enumerate(signal_names):
        signal = mdf.get(signal_name)

        ax = rows[i][0]
        ax.set_title(signal_name)
        ax.plot(signal.timestamps, signal.samples)

        signal_name_sanitized = signal_name.replace('Model Root/', '').replace('/In1', '')
        signal_map[signal_name_sanitized] = signal

    e.commit_fig('all.pdf', fig)

    # - EXTRACTING THE DIFFERENT MEASUREMENTS

    e.info('extracting the different measurements from the file by looking for switch == LOW')
    timestamps = signal.timestamps
    measurements = []
    measurement_data = defaultdict(list)
    start_time = None
    for i in range(len(timestamps)):
        t = timestamps[i]

        if start_time is None and signal_map['switch_mess'][i] < 0.1:
            start_time = t

        if start_time is not None and signal_map['switch_mess'][i] < 0.1:
            measurement_data['timestamps'].append(t - start_time)
            for name, signal in signal_map.items():
                if name == 'y_out_mess':
                    value = 1.3 - signal.samples[i]
                else:
                    value = signal.samples[i]

                measurement_data[name].append(value)

        if start_time is not None and signal_map['switch_mess'][i] > 0.9:
            measurements.append(dict(measurement_data))
            measurement_data = defaultdict(list)
            start_time = None

    e.info(f'extracted a total of {len(measurements)} measurements from the file')

    e.info('plotting the different measurements...')
    pdf_path = os.path.join(e.path, 'measurements.pdf')
    with PdfPages(pdf_path) as pdf:

        for m, measurement_data in enumerate(measurements):

            n_rows = len(measurement_data)
            fig, rows = plt.subplots(ncols=1, nrows=n_rows, squeeze=False, figsize=(16, 4 * n_rows))
            fig.suptitle(f'measurement {m}')

            signal_map = {}
            for i, (key, values) in enumerate(measurement_data.items()):
                ax = rows[i][0]
                ax.set_title(key)
                ax.plot(measurement_data['timestamps'], measurement_data[key])

            pdf.savefig(fig)
            plt.close(fig)

    e.commit_json('measurements.json', measurements)
