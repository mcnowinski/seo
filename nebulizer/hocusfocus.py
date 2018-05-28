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

# for plots
import matplotlib
matplotlib.use('Agg')  # don't need display
import matplotlib.pyplot as plt
from astropy.visualization import astropy_mpl_style
plt.style.use(astropy_mpl_style)
plt.ioff()

# user, hardcode for now
user = 'hocusfocus'

#
# main
#

# set up logger
logger = log.get_logger('hocusfocus')

# simulate? set to True
simulate = False

# time between checks for object observability in seconds
delay_time = 30

# image exposure in seconds
exposure = 30

# image binning
binning = 2

# image filter
filter = 'clear'

# minimum altitude to observe this target
min_obs_alt = 60.0

# list of references stars for focus tests
ref_stars_fname = '/home/mcnowinski/seo/nebulizer/SDSS_Standard_Stars'

# path to sextractor executable, configuration, and output files
sextractor_bin_path = '/home/mcnowinski/sex/sextractor/bin/sex'
sextractor_sex_path = '/home/mcnowinski/seo/nebulizer/hf.sex'
sextractor_param_path = '/home/mcnowinski/seo/nebulizer/hf.param'
sextractor_conv_path = '/home/mcnowinski/seo/nebulizer/hf.conv'
sextractor_cat_path = '/home/mcnowinski/seo/nebulizer/hf.cat'

# path to psfex executable, configuration, and output file
psfex_bin_path = '/home/mcnowinski/psfex/bin/psfex'
psfex_cfg_path = '/home/mcnowinski/seo/nebulizer/hf.psfex'
psfex_psf_path = '/home/mcnowinski/seo/nebulizer/hf.psf'

# focus plot
plt_path = '/tmp/hocusfocus.png'

# image path
image_path = '/tmp'

# image filename
image_filename = 'hocusfocus.fits'

# seo
observatory = Observatory('G52', 38.2886, -122.50400, 8, 'US/Pacific')

# seo telescope
telescope = Telescope(simulate)

# list of observations
observations = []

# list of image stacks
stacks = []

# initialize array covering a range of focus positions
pass1_array = [4650, 4675, 4700, 4725, 4750, 4775, 4800,
               4825, 4850, 4875, 4900, 4925, 4950, 4975, 5000]
pass1_array_focus = np.zeros((len(pass1_array), 2))

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

# wait for sun to set
telescope.checkSun(True)

# find the best target (currently highest in the sky)
# this is overly complicated, but I know it works ;)
# start up the scheduler
scheduler = Scheduler(observatory, observations)
telescope.slackdebug('Identifying calibration star...')
target_observation = scheduler.whatsHighest()

if target_observation:
    logger.debug(target_observation.toString())
else:
    logger.error('No target found. Aborted.')
    sys.exit(1)

telescope.slackdebug('Calibration star is %s.' %
                     target_observation.target.getName())

# check clouds
telescope.checkClouds()

# open up the telescope
telescope.crackit()

# point the scope
telescope.slackdebug('Pointing telescope to %s...' %
                     target_observation.target.getName())
telescope.pinpointier(target_observation)

# get current focus position in case something goes awry
telescope.slackdebug('Original focus position is %s.' % telescope.getFocus())
focus_position_default = int(telescope.getFocus())

# calculate PSF for pass1_array focus positions
for i, focus_position in enumerate(pass1_array):
    # set focus to ith position of pass1_array
    telescope.slackdebug('Setting focus position to %s...' % focus_position)
    telescope.setFocus(focus_position)

    # take image
    telescope.slackdebug('Taking calibration image...')
    telescope.image(exposure, filter, binning, image_path, image_filename)
    telescope.slackpreview(image_path+'/'+image_filename)

    # perform photometry
    telescope.slackdebug('Performing photometry on image...')
    (output, error, pid) = telescope.runSubprocess(
        [sextractor_bin_path, image_path+'/'+image_filename, '-c', sextractor_sex_path])

    # perform psf extraction
    telescope.slackdebug('Estimating point spread function of image...')
    (output, error, pid) = telescope.runSubprocess(
        [psfex_bin_path, sextractor_cat_path, '-c', psfex_cfg_path])

    psf = fits.open(psfex_psf_path)
    fwhm = psf[1].header['PSF_FWHM']

    pass1_array_focus[i] = focus_position, fwhm

    telescope.slackdebug(
        'For a focus position of %s, estimated FWHM is %s.' % (focus_position, fwhm))
    logger.info('focus_position=%s, fwhm=%s' % (focus_position, fwhm))

    telescope.pinpointier(target_observation, False)

# fit data points to a 2nd-deg polynomial
pass1_fit = np.polyfit(pass1_array_focus[:, 0], pass1_array_focus[:, 1], 2)
pass1_fit_focus_pos = int(-pass1_fit[1]/(2*pass1_fit[0]))

# save focus pos array
#np.savetxt("/home/sirius/focus/"+folder+"/"+folder+".dat", pass1_array_focus, fmt='%.5f', header='focus_pos PSF')

# graph focus fits
array = np.array(pass1_array_focus)
plt.scatter(array[:, 0], array[:, 1])
x = np.arange(4000, 5100)
y = pass1_fit[0]*x**2+pass1_fit[1]*x+pass1_fit[2]
plt.ylim(2, 5.5)
plt.xlim(4550, 5100)
plt.xlabel('Focus Position')
plt.ylabel('FWHM')
plt.savefig(plt_path, bbox_inches='tight')
plt.close()

telescope.slackimage(plt_path)

# set focus to minimum
telescope.slackdebug('Setting final focus position to %d...' %
                     pass1_fit_focus_pos)
telescope.setFocus(pass1_fit_focus_pos)

telescope.squeezeit()

telescope.slackdebug('Hocusfocus complete.')
