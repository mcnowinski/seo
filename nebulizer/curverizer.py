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
user = 'curverizer'

#
# main
#

# set up logger
logger = log.get_logger('curverizer')

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

input = ['mcnowinski', 'Solar System', '3982', 120, 1, 'clear', 2, Sequence.CONTINUOUS]

user = input[0].replace('@', '.')
obstype = input[1]
name = input[2]
exposure = input[3]
count = input[4]
filters = input[5].split(',')
binning = input[6]
repeat = input[7];
if obstype == "Solar System":
    objects = scheduler.findSolarSystemObjects(name.strip().upper())
else:
    objects = scheduler.findObjects(name.strip().upper())
if len(objects) == 0:
    logger.error('Could not find matching object for %s.' % name)
    sys.exit(1)
else:
    if len(objects) > 1:
        logger.error('Found multiple matching objects for %s.' % name)
        logger.error('Will take the first in the list (%s).' % name)
        sys.exit(1)
    target = Target(objects[0]['name'], objects[0]
                    ['RA'], objects[0]['DEC'])
    stacks = []
    for filter in filters:
        filter = filter.strip()
        # skip the darks for now
        if filter.lower() != 'dark':
            stacks.append(Stack(float(exposure), filter,
                                int(binning), int(count)))
    observations.append(Observation(observatory, target,
                                    Sequence(stacks, repeat), min_obs_alt, user))

telescope.slackdebug('Starting curverizer...')

# add observations
scheduler.addObservations(observations)

# loop through the observations
asteroid_observation = scheduler.whatsInBetween()

if asteroid_observation:
    logger.debug(asteroid_observation.toString())
else:
    logger.error('Target (%s) is not observable.' %
                 asteroid_observation.target.getName())
    sys.exit(1)

# wait for sun to set
telescope.checkSun(True)

# wait for target to be available
while Time.now() < asteroid_observation.min_obs_time:
    wait_time_s = (asteroid_observation.min_obs_time-Time.now()).sec
    telescope.slackdebug('The observation (%s) will start in %d min (at %s)...' % (
        asteroid_observation.target.getName(), wait_time_s/60, asteroid_observation.min_obs_time.iso[:-7]))

# start observations
telescope.crackit()
telescope.slackdebug("Observations for %s are starting..." %
                     (asteroid_observation.target.getName()))
# point the scope
telescope.pinpointier(asteroid_observation)

if asteroid_observation.max_alt_time > asteroid_observation.min_obs_time:
    # get images of target until time of max altitude time
    while Time.now() < asteroid_observation.max_alt_time:
        # check sun, clouds, slit, etc.
        telescope.checkSun()
        if telescope.is_cracked:
            telescope.checkSlit()
        telescope.checkAlt()
        telescope.checkClouds()

        telescope.getImageStacks(asteroid_observation, None, False)

# get images of target until target has set
while Time.now() < asteroid_observation.max_obs_time:
    # check sun, clouds, slit, etc.
    telescope.checkSun()
    if telescope.is_cracked:
        telescope.checkSlit()
    telescope.checkAlt()
    telescope.checkClouds()

    telescope.getImageStacks(asteroid_observation, None, False)

telescope.slackdebug("Observations for %s are complete." %
                     (asteroid_observation.target.getName()))

telescope.done()
