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

# user, hardcode for now
user = 'observerizer'

#
# main
#

# set up logger
logger = log.setup_custom_logger('uc')

# list of observations
input_fname = '/home/mcnowinski/seo/nebulizer/observations.txt'

# simulate? set to True
simulate = False

# time between checks for object observability in seconds
delay_time = 30

# min time available for background observations in seconds
min_background_time = 10*60  # 10 minutes

# seo
observatory = Observatory('G52', 38.2886, -122.50400, 8, 'US/Pacific')

# seo telescope
telescope = Telescope(simulate)

# default min observation alt in degrees
min_obs_alt = 30.0

# list of observations
observations = []

# init the scheduler
scheduler = Scheduler(observatory)

# read in observations from file
with open(input_fname) as f:
    inputs = f.readlines()
inputs = [x.strip() for x in inputs if not x.startswith('#')]

# jmartinez8@uchicago.edu	General	NGC 3941	240	1	clear, dark, u-band, g-band, r-band, i-band	2	No
num_targets = 0
for input in inputs:
    input = input.split('\t')
    user = input[0].replace('@', '.')
    name = input[2]
    exposure = input[3]
    count = input[4]
    filters = input[5].split(',')
    binning = input[6]
    objects = scheduler.findObjects(name.strip().upper())
    if len(objects) == 0:
        print 'Error. Could not find matching object for %s.' % name
    elif len(objects) > 1:
        print 'Error. Found multiple matching objects for %s.' % name
    else:
        num_targets += 1
        target = Target(objects[0]['name'], objects[0]
                        ['RA'], objects[0]['DEC'])
        stacks = []
        for filter in filters:
            filter = filter.strip()
            #skip the darks for now
            if filter.lower() == 'dark':
                stacks.append(Stack(float(exposure), filter,
                                    int(binning), int(count), False))
        observations.append(Observation(observatory, target,
                                        Sequence(stacks, 1), min_obs_alt, user))

print 'Successfully identified %d of %d targets.' % (
    num_targets, len(observations))

telescope.slackdebug('Starting dark observerizer...')

for observation in observations:
    telescope.getImageStacks(observation, None, False)

telescope.slackdebug(
    'Dark Observerization complete. Performed %d observations.' % count)
