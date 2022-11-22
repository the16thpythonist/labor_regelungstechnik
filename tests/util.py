import os
import sys
import logging
import pathlib

PATH = pathlib.Path(__file__).parent.absolute()
ASSETS_PATH = os.path.join(PATH, 'assets')

LOG = logging.Logger("TESTING")
LOG.addHandler(logging.StreamHandler(sys.stdout))
