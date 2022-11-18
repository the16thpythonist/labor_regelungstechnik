import os
import re
from pprint import pprint, PrettyPrinter

import numpy as np
import matplotlib.pyplot as plt
from pycomex.experiment import Experiment
from pycomex.util import Skippable
from asammdf import MDF

from labor_regelungstechnik.utils import MEASUREMENTS_PATH

MEASUREMENTS_FILE_NAME = 'messungen_17_11_22.mf4'

BASE_PATH = os.getcwd()
NAMESPACE = 'extract_measurements'
DEBUG = True
with Skippable(), (e := Experiment(base_path=BASE_PATH, namespace=NAMESPACE, glob=globals())):

    measurements_file_path = os.path.join(MEASUREMENTS_PATH, MEASUREMENTS_FILE_NAME)
    e.info(f'starting to extract measurements from recorded file "{measurements_file_path}"')

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

    for i, signal_name in enumerate(signal_names):
        signal = mdf.get(signal_name)

        ax = rows[i][0]
        ax.set_title(signal_name)
        ax.plot(signal.timestamps, signal.samples)

    e.commit_fig('plots.pdf', fig)
