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

# path to configuration
cfg_path = "/home/mcnowinski/seo/nebulizer/curvacious.json"

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
wait_time_s = cfg['wait_time_s']

# build main asteroid observation
observation_json = cfg['observations']['asteroid']
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
        stack_json['binning']), int(stack_json['count']))
    logger.debug(stack.toString().replace('\n', '; '))
    stacks.append(stack)
# build sequence
sequence = Sequence(stacks, int(sequence_json['repeat']))
logger.debug(sequence.toString().replace('\n', '; '))
# build main observation
asteroid_main_observation = Observation(
    observatory, target, sequence, min_obs_alt, user)
# get min, max, and max alt obs times
asteroid_main_observation.getTimes()
logger.debug(asteroid_main_observation.toString().replace('\n', '; '))

# build bracket asteroid observation
sequence_json = observation_json['sequences']['bracket']
stacks_json = sequence_json['stacks']
# build image stacks
stacks = []
for stack_json in stacks_json:
    stack = Stack(float(stack_json['exposure']), stack_json['filters'], int(
        stack_json['binning']), int(stack_json['count']))
    logger.debug(stack.toString().replace('\n', '; '))
    stacks.append(stack)
# build sequence
sequence = Sequence(stacks, int(sequence_json['repeat']))
logger.debug(sequence.toString().replace('\n', '; '))
# build bracket observation
asteroid_bracket_observation = Observation(
    observatory, target, sequence, min_obs_alt, user)
logger.debug(asteroid_bracket_observation.toString().replace('\n', '; '))

# build referece (star) observation
observation_json = cfg['observations']['reference']
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
        stack_json['binning']), int(stack_json['count']))
    logger.debug(stack.toString().replace('\n', '; '))
    stacks.append(stack)
# build sequence
sequence = Sequence(stacks, int(sequence_json['repeat']))
logger.debug(sequence.toString().replace('\n', '; '))
# build bracket observation
reference_observation = Observation(
    observatory, target, sequence, min_obs_alt, user)
# get min, max, and max alt obs times
reference_observation.getTimes()    
logger.debug(reference_observation.toString().replace('\n', '; '))

# start observations
telescope.slackdebug('Starting curvacious...')

# wait for sun to set
telescope.checkSun(True)

# wait for target to be available
while Time.now() < asteroid_main_observation.min_obs_time:
    wait_time_s = (asteroid_main_observation.min_obs_time-Time.now()).sec
    telescope.slackdebug('The observation (%s) will start in %d min (at %s)...' % (
        asteroid_main_observation.target.getName(), wait_time_s/60, asteroid_main_observation.min_obs_time.iso[:-7]))
    time.sleep(wait_time_s)

# start observations
telescope.slackdebug("Starting observations...")
telescope.crackit()

# do asteroid observations while (if) we need to wait for reference to be availabe
if Time.now() < reference_observation.min_obs_time:
    telescope.slackdebug("Performing asteroid (%s) observations until %s..." % (
        asteroid_main_observation.target.getName(), reference_observation.min_obs_time.iso[:-7]))
while Time.now() < reference_observation.min_obs_time:
    # check sun, clouds, slit, etc.
    telescope.checkSun()
    if telescope.is_cracked:
        telescope.checkSlit()
    telescope.checkAlt()
    telescope.checkClouds()
    telescope.getImageStacks(asteroid_main_observation, None, False)

# do reference observations
telescope.slackdebug("Performing reference observations...")
telescope.slackdebug("Observing %s..." %
                     (reference_observation.target.getName()))
# check sun, clouds, slit, etc.
telescope.checkSun()
if telescope.is_cracked:
    telescope.checkSlit()
telescope.checkAlt()
telescope.checkClouds()
# point the scope
telescope.pinpointier(reference_observation)
telescope.getImage(reference_observation)
telescope.slackdebug("Reference observations complete.")

# do reference observations
telescope.slackdebug("Performing bracket observations...")
telescope.slackdebug("Observing %s..." %
                     (asteroid_bracket_observation.target.getName()))
# check sun, clouds, slit, etc.
telescope.checkSun()
if telescope.is_cracked:
    telescope.checkSlit()
telescope.checkAlt()
telescope.checkClouds()
# point the scope
telescope.pinpointier(asteroid_bracket_observation)
telescope.getImage(asteroid_bracket_observation)
telescope.slackdebug("Bracket observations complete.")

telescope.slackdebug("Performing asteroid (%s) observations..." %
                     (asteroid_main_observation.target.getName()))
# get images of target until target has set
while Time.now() < asteroid_main_observation.max_obs_time:
    # check sun, clouds, slit, etc.
    telescope.checkSun()
    if telescope.is_cracked:
        telescope.checkSlit()
    telescope.checkAlt()
    telescope.checkClouds()
    telescope.getImageStacks(asteroid_main_observation, None, False)

telescope.slackdebug("All observations for %s are complete." %
                     (asteroid_main_observation.target.getName()))

telescope.done()
