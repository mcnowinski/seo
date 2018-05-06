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
user = 'asterizer'

#
# main
#

# set up logger
logger = log.setup_custom_logger('asterizer')

# list of observations
input_fname = '/home/mcnowinski/seo/nebulizer/asteroids.txt'

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

# user,object,exp,bin,filter,count
num_targets = 0
for input in inputs:
    input = input.split(',')
    user = input[0]
    name = input[1]
    exposure = input[2]
    binning = input[3]
    filters = input[4].split(' ')
    count = int(input[5])
    objects = scheduler.findSolarSystemObjects(name.strip().upper())
    if len(objects) == 0:
        print 'Error. Could not find matching object for %s.' % name
    elif len(objects) > 1:
        print 'Error. Found multiple matching objects for %s.' % name
    else:
        num_targets += 1
        target = Target(objects[0]['name'], objects[0]
                        ['RA'], objects[0]['DEC'])
        stacks = []
        for i in range(0, count):
            for filter in filters:
                filter = filter.strip()
                # skip the darks for now
                if filter.lower() != 'dark':
                    stacks.append(Stack(float(exposure), filter,
                                        int(binning), 1))
        observations.append(Observation(observatory, target,
                                        Sequence(stacks, 1), min_obs_alt, user))

print 'Successfully identified %d of %d targets.' % (
    num_targets, len(inputs))

telescope.slackdebug('Starting observerizer...')

# add observations
scheduler.addObservations(observations)

# loop through the observations
next_observation = scheduler.whatsNext()

# for obs in observations:
#    print obs.toString()

if next_observation:
    print next_observation.toString()

# wait for sun to set
telescope.checkSun(True)

count = 0
while next_observation != None:

    # check sun, clouds, slit, etc.
    telescope.checkSun()
    if telescope.is_cracked:
        telescope.checkSlit()
    telescope.checkAlt()
    telescope.checkClouds()

    if Time.now() >= next_observation.min_obs_time:  # its time to look at this object
        count += 1
        if not telescope.is_cracked:
            telescope.crackit()
        # perform observation
        telescope.slackdebug("Observations for %s are starting..." %
                             (next_observation.target.getName()))
        # point the scope
        telescope.pinpointier(next_observation)
        # get the images
        telescope.getImage(next_observation)
        telescope.slackdebug("Observations for %s are complete." %
                             (next_observation.target.getName()))
        scheduler.isDone(next_observation)  # mark as complete
    else:
        wait_time_s = (next_observation.min_obs_time-Time.now()).sec
        telescope.slackdebug('The next observation (%s) will start in %d min (at %s)...' % (
            next_observation.target.getName(), wait_time_s/60, next_observation.min_obs_time.iso[:-7]))
        time.sleep(delay_time)

    # if its been a while, check our target list again
    next_observation = scheduler.whatsNext()

telescope.slackdebug(
    'Asterizer complete. Performed %d observations.' % count)
