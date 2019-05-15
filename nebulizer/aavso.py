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
cfg_path = "/home/mcnowinski/seo/nebulizer/aavso.json"

# set up logger
logger = log.get_logger('aavso')

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

# build main star observation
observation_json = cfg['observations']['star']
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
star_main_observation = Observation(
    observatory, target, sequence, min_obs_alt, user)
# get min, max, and max alt obs times
star_main_observation.getTimes()
logger.debug(star_main_observation.toString().replace('\n', '; '))

# start observations
telescope.slackdebug('Starting aavso...')

# wait for sun to set
telescope.checkSun(True)

# wait for target to be available
while Time.now() < star_main_observation.min_obs_time:
    time_left_s = (star_main_observation.min_obs_time-Time.now()).sec
    telescope.slackdebug('The observation (%s) will start in %d min (at %s)...' % (
        star_main_observation.target.getName(), time_left_s/60, star_main_observation.min_obs_time.iso[:-7]))
    time.sleep(wait_time_s)

# start observations
telescope.slackdebug("Starting observations...")
telescope.crackit()

telescope.slackdebug("Performing star (%s) observations..." %
                     (star_main_observation.target.getName()))
# get images of target until target has set
while Time.now() < star_main_observation.max_obs_time:
    # check sun, clouds, slit, etc.
    telescope.checkSun()
    if telescope.is_cracked:
        telescope.checkSlit()
    telescope.checkAlt()
    telescope.checkClouds()
    telescope.getImageStacks(star_main_observation, None, False)

telescope.slackdebug("Observations for %s are complete." %
                     (star_main_observation.target.getName()))

telescope.done()
