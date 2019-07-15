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


def doObservations():
    telescope.slackdebug("Performing observations for %s..." %
                         (target_main_observation.target.getName()))
    # check sun, clouds, slit, etc.
    telescope.checkSun()
    if telescope.is_cracked:
        telescope.checkSlit()
    telescope.checkAlt()
    telescope.checkClouds()
    # point the scope
    telescope.pinpoint(target_main_observation)
    telescope.getImage(target_main_observation)
    telescope.slackdebug("Observations complete.")


# path to configuration
cfg_path = "/home/mcnowinski/seo/nebulizer/transformation.json"

# set up logger
logger = log.get_logger('transformerizer')

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

# build observation
observation_json = cfg['observations']
target_json = observation_json['target']
sequence_json = observation_json['sequences']
stacks_json = sequence_json['stacks']
# build target
target = Target.from_name(
    target_json['name'], observatory, target_json['type'], target_json.get('ra_offset'), target_json.get('dec_offset'))
logger.debug(target.toString().replace('\n', '; '))
# build image stacks
stacks = []
for stack_json in stacks_json:
    stack = Stack(float(stack_json['exposure']), stack_json['filters'], int(
        stack_json['binning']), int(stack_json['count']), stack_json['do_pinpoint'] if 'do_pinpoint' in stack_json else True)
    logger.debug(stack.toString().replace('\n', '; '))
    stacks.append(stack)
# build sequence
sequence = Sequence(stacks, int(
    sequence_json['repeat']), sequence_json['do_pinpoint'] if 'do_pinpoint' in sequence_json else True)
logger.debug(sequence.toString().replace('\n', '; '))
# build main observations
target_main_observation = Observation(
    observatory, target, sequence, min_obs_alt, user)
# get min, max, and max alt obs times
target_main_observation.getTimes()
logger.debug(target_main_observation.toString().replace('\n', '; '))

# start observations
telescope.slackdebug('Starting transformerizer...')

# wait for sun to set
telescope.checkSun(True)

# wait for target to be available
while Time.now() < target_main_observation.min_obs_time:
    wait_time_s = (target_main_observation.min_obs_time-Time.now()).sec
    telescope.slackdebug('The observation (%s) will start in %d min (at %s)...' % (
        target_main_observation.target.getName(), wait_time_s/60, target_main_observation.min_obs_time.iso[:-7]))
    time.sleep(delay_time)

# start observations
telescope.crackit()

doObservations()

telescope.done()
