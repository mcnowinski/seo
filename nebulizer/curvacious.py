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
observatory = Observatory(cfg['observatory']['code'], cfg['observatory']['latitude'], cfg['observatory']['longitude'], cfg['observatory']['altitude'], cfg['observatory']['timezone'])

# init seo telescope
telescope = Telescope(cfg['simulate'])

# pause time while waiting for object to become availab.e
wait_time_s = cfg['wait_time_s']

# build main observation
observation_json = cfg['observations']['main']
target_json = cfg['observations']['main']['target']
sequence_json = cfg['observations']['main']['sequence']
stacks_json = cfg['observations']['main']['sequence']['stacks']
# build target
target = Target.from_name(target_json['name'], observatory, target_json['type'])
logger.debug(target.toString().replace('\n','; '))
# build image stacks
stacks = []
for stack_json in stacks_json:
    stack = Stack(float(stack_json['exposure']), stack_json['filters'], int(stack_json['binning']), int(stack_json['count']))
    logger.debug(stack.toString().replace('\n','; '))
    stacks.append(stack)
# build sequence
sequence = Sequence(stacks, int(sequence_json['repeat']))
logger.debug(sequence.toString().replace('\n','; '))
# build main observation
main_observation = Observation(observatory, target, sequence, min_obs_alt, user)
main_observation.getTimes()
logger.debug(main_observation.toString().replace('\n','; '))

#build comparison observations
comparison_observations = []
observations_json = cfg['observations']['comparisons']
for observation_json in observations_json:
    target_json = observation_json['target']
    sequence_json = observation_json['sequence']
    stacks_json = observation_json['sequence']['stacks']  
    # build target  
    target = Target.from_name(target_json['name'], observatory, target_json['type'])
    logger.debug(target.toString().replace('\n','; '))
    # build image stacks
    stacks = []
    for stack_json in stacks_json:
        stack = Stack(float(stack_json['exposure']), stack_json['filters'], int(stack_json['binning']), int(stack_json['count']))
        logger.debug(stack.toString().replace('\n','; '))
        stacks.append(stack)
    # build sequence
    sequence = Sequence(stacks, int(sequence_json['repeat']))
    logger.debug(sequence.toString().replace('\n','; '))
    # build main observation
    comparison_observation = Observation(observatory, target, sequence, min_obs_alt, user)
    comparison_observation.getTimes()
    logger.debug(comparison_observation.toString().replace('\n','; '))
    comparison_observations.append(comparison_observation)

#start observations
telescope.slackdebug('Starting curvacious...')

# wait for sun to set
telescope.checkSun(True)

# wait for target to be available
while Time.now() < main_observation.min_obs_time:
    wait_time_s = (main_observation.min_obs_time-Time.now()).sec
    telescope.slackdebug('The observation (%s) will start in %d min (at %s)...' % (
        main_observation.target.getName(), wait_time_s/60, main_observation.min_obs_time.iso[:-7]))
    time.sleep(wait_time_s)

# start observations
telescope.crackit()
telescope.slackdebug("Observations are starting...")

# do comparison observations
total_time_comparison_observations_s = 0
if len(comparison_observations) > 0:
    total_time_comparison_observations_s = 15*60
    telescope.slackdebug("Starting comparison observations...")
    for comparison_observation in comparison_observations:
        telescope.slackdebug("Observing %s..." %
                            (comparison_observation.target.getName()))
         # check sun, clouds, slit, etc.
        telescope.checkSun()
        if telescope.is_cracked:
            telescope.checkSlit()
        telescope.checkAlt()
        telescope.checkClouds()       
        # point the scope
        telescope.pinpointier(comparison_observation)
        for repeat in range(0, comparison_observation.sequence.repeat):
            telescope.getImageStacks(comparison_observation)
    telescope.slackdebug("Comparison observations complete.")    

telescope.slackdebug("Starting target observations...")
telescope.slackdebug("Observing %s..." %
                    (main_observation.target.getName()))
# get images of target until target has set
while Time.now() + TimeDelta(total_time_comparison_observations_s, format='sec') < main_observation.max_obs_time:
    # check sun, clouds, slit, etc.
    telescope.checkSun()
    if telescope.is_cracked:
        telescope.checkSlit()
    telescope.checkAlt()
    telescope.checkClouds()
    telescope.getImageStacks(main_observation, None, False)
telescope.slackdebug("Target observations complete.") 

if len(comparison_observations) > 0:
    telescope.slackdebug("Starting comparison observations...")
    for comparison_observation in comparison_observations:
        telescope.slackdebug("Observing %s..." %
                            (comparison_observation.target.getName()))
         # check sun, clouds, slit, etc.
        telescope.checkSun()
        if telescope.is_cracked:
            telescope.checkSlit()
        telescope.checkAlt()
        telescope.checkClouds()       
        # point the scope
        telescope.pinpointier(comparison_observation)
        for repeat in range(0, comparison_observation.sequence.repeat):
            telescope.getImageStacks(comparison_observation)
    telescope.slackdebug("Comparison observations complete.")  

telescope.slackdebug("all observations for %s are complete." %
                     (asteroid_observation.target.getName()))

telescope.done()







