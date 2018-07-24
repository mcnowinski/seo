import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)

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
from astropy.time import TimeDelta
import time
import json


def calibrationObservations():
    telescope.slackdebug("Performing calibration observations for %s..." %
                         (asteroid_main_observation.target.getName()))
    # check sun, clouds, slit, etc.
    telescope.checkSun()
    if telescope.is_cracked:
        telescope.checkSlit()
    telescope.checkAlt()
    telescope.checkClouds()
    # point the scope
    telescope.pinpointier(asteroid_calibration_observation)
    telescope.getImage(asteroid_calibration_observation)
    telescope.slackdebug("Calibration observations complete.")


# path to configuration
cfg_path = "/home/mcnowinski/seo/nebulizer/curvacious2.json"

# set up logger
logger = log.get_logger('curvacious')

# load target and comparison observations
with open(cfg_path) as f:
    cfg = json.load(f)

# user, hardcode for now
user = cfg['user']

# min obs altitude
min_obs_alt = float(cfg['min_obs_alt'])

# seo
observatory = Observatory(cfg['observatory']['code'], cfg['observatory']['latitude'], cfg['observatory']
                          ['longitude'], cfg['observatory']['altitude'], cfg['observatory']['timezone'])

# init seo telescope
telescope = Telescope(cfg['simulate'])

# pause time while waiting for object to become available
delay_time = cfg['delay_time']

# build main asteroid observation
observation_json = cfg['observations']
target_json = observation_json['target']
sequence_json = observation_json['sequences']['main']
stacks_json = sequence_json['stacks']
# build target
target = Target.from_name(
    target_json['name'], observatory, target_json['type'])
logger.debug(target.toString().replace('\n', '; '))
# build image stacks
stacks = []
for stack_json in stacks_json:
    stack = Stack(float(stack_json['exposure']), stack_json['filters'], int(
        stack_json['binning']), int(stack_json['count']), stack_json['do_pinpoint'] if 'do_pinpoint' in stack_json else True)
    logger.debug(stack.toString().replace('\n', '; '))
    stacks.append(stack)
# build sequence
sequence = Sequence(stacks, int(sequence_json['repeat']))
logger.debug(sequence.toString().replace('\n', '; '))
# build main observations
asteroid_main_observation = Observation(
    observatory, target, sequence, min_obs_alt, user)
# get min, max, and max alt obs times
asteroid_main_observation.getTimes()
logger.debug(asteroid_main_observation.toString().replace('\n', '; '))

# build calibration asteroid/star observations
sequence_json = observation_json['sequences']['calibration']
stacks_json = sequence_json['stacks']
# build image stacks
stacks = []
for stack_json in stacks_json:
    stack = Stack(float(stack_json['exposure']), stack_json['filters'], int(
        stack_json['binning']), int(stack_json['count']), stack_json['do_pinpoint'] if 'do_pinpoint' in stack_json else True)
    logger.debug(stack.toString().replace('\n', '; '))
    stacks.append(stack)
# build sequence
sequence = Sequence(stacks, int(sequence_json['repeat']))
logger.debug(sequence.toString().replace('\n', '; '))
# build calibration observations
asteroid_calibration_observation = Observation(
    observatory, target, sequence, min_obs_alt, user)
asteroid_calibration_observation_duration_s = sequence.getDuration()
logger.debug(asteroid_calibration_observation.toString().replace('\n', '; '))

# start observations
telescope.slackdebug('Starting curvacious...')

# wait for sun to set
telescope.checkSun(True)

# wait for target to be available
while Time.now() < asteroid_main_observation.min_obs_time:
    wait_time_s = (asteroid_main_observation.min_obs_time-Time.now()).sec
    telescope.slackdebug('The observation (%s) will start in %d min (at %s)...' % (
        asteroid_main_observation.target.getName(), wait_time_s/60, asteroid_main_observation.min_obs_time.iso[:-7]))
    time.sleep(delay_time)

# start observations
telescope.slackdebug("Starting observations...")
telescope.crackit()

calibrationObservations()

# get images of target until target reaches max altitude, if it hasn't already
doCalibrationObservations = False
while Time.now() + TimeDelta(asteroid_calibration_observation_duration_s, format='sec') < asteroid_main_observation.max_alt_time:
    telescope.slackdebug("Performing main observations for %s..." %
                         (asteroid_main_observation.target.getName()))
    # check sun, clouds, slit, etc.
    telescope.checkSun()
    if telescope.is_cracked:
        telescope.checkSlit()
    telescope.checkAlt()
    telescope.checkClouds()
    telescope.getImageStacks(asteroid_main_observation, None, False)
    doCalibrationObservations = True

if doCalibrationObservations:
    calibrationObservations()

# get images of target until target has set
doCalibrationObservations = False
while Time.now() + TimeDelta(asteroid_calibration_observation_duration_s, format='sec') < asteroid_main_observation.max_obs_time:
    telescope.slackdebug("Performing main observations for %s..." %
                         (asteroid_main_observation.target.getName()))
    # check sun, clouds, slit, etc.
    telescope.checkSun()
    if telescope.is_cracked:
        telescope.checkSlit()
    telescope.checkAlt()
    telescope.checkClouds()
    telescope.getImageStacks(asteroid_main_observation, None, False)
    doCalibrationObservations = True

if doCalibrationObservations:
    calibrationObservations()

telescope.slackdebug("All observations for %s are complete." %
                     (asteroid_main_observation.target.getName()))

telescope.done()
