# import classes from the chultun module
from chultun import Target  # name, ra, dec
from chultun import Stack  # exposure, filter, binning, count
from chultun import Sequence  # stacks, repeat
from chultun import Observatory  # code, latitude, longitude, altitude
from chultun import Observation  # target, sequence
from chultun import Scheduler  # observatory, observations
from chultun import Telescope  # the telescope commands

import log
import datetime
import traceback
import sys
import subprocess
from astropy.time import Time
import time
import numpy as np
from astropy.io import fits

# seo telescope
telescope = Telescope(False)

telescope.setFilter('sii')
