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
user = 'chutlun'

#
# main
#

# set up logger
logger = log.setup_custom_logger('hocusfocus')

# simulate? set to True
simulate = True

# time between checks for object observability in seconds
delay_time = 30

# image exposure in seconds
exposure = 30

# image binning
binning = 2

# image filter
filter = 'clear'

# minimum altitude to observe this target
min_obs_alt = 30.0

# min time available for background observations in seconds
min_background_time = 10*60  # 10 minutes

# list of references stars for focus tests
ref_stars_fname = '/home/mcnowinski/seo/nebulizer/SDSS_Standard_Stars'

# seo
observatory = Observatory('G52', 38.2886, -122.50400, 8, 'US/Pacific')

# seo telescope
telescope = Telescope(simulate)

# list of observations
observations = []

# list of image stacks
stacks = []

# read in reference stars from file
with open(ref_stars_fname) as f:
    ref_stars = f.readlines()
ref_stars = [x.strip() for x in ref_stars if not x.startswith('#')]

# create a list of observations from the reference stars
stack = Stack(exposure, filter, binning, 1)
stacks = [stack]
sequence = Sequence(stacks, 1)
for ref_star in ref_stars:
    ref_star_data = ref_star.split()
    name = ref_star_data[0]
    ra = ref_star_data[1]
    dec = ref_star_data[2]
    target = Target(name, ra, dec)
    observations.append(Observation(
        observatory, target, sequence, min_obs_alt, user))

telescope.slackdebug('Starting hocusfocus...')

# start up the scheduler
scheduler = Scheduler(observatory, observations)

# loop through the observations
next_observation = scheduler.whatsNext()

if next_observation:
    print next_observation.toString()

telescope.slackdebug('Hocusfocus complete.')
