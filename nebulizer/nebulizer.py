#import classes from the chultun module
from chultun import Target #name, ra, dec
from chultun import Stack #exposure, filter, binning, count
from chultun import Sequence #stacks, repeat
from chultun import Observatory #code, latitude, longitude, altitude
from chultun import Observation #target, sequence
from chultun import Scheduler #observatory, observations
from chultun import Telescope #the telescope commands

import log
import datetime
import traceback
import sys
import subprocess
from astropy.time import Time
import time

#user, hardcode for now
user='chutlun'
        
#
#main
#

#set up logger
logger = log.setup_custom_logger('nebulizer')

#simulate? set to True
simulate = False

#time between checks for object observability in seconds
delay_time = 30

#seo
observatory = Observatory('G52', 38.2886, -122.50400, 8, 'US/Pacific')

#seo telescope
telescope = Telescope(simulate)

#list of observations
observations = []
#list of image stacks
stacks = []

#NAME,      RA (ICRS),      DEC (ICRS),     B FLUX, V FLUX, R FLUX, J FLUX, H FLUX, K FLUX, NOTES
#PN M 1-4   03 41 43.439    +52 16 59.85    15.86,14.27,15.25,12.529,12.204,11.245,
#PN K 3-63  21 39 11.976    +55 46 03.94    ,,,,,,I did 120 s for this one and it seemed okay
#PN K 3-66  04 36 37.23     +33 39 30.0     14.8,,12.5,13.205,12.957,12.353,
#IC 2149    05 56 23.862    +46 06 17.50    10.5,10.2,,10.372,10.308,9.699,
#PN M 1-6   06 35 45.126    -00 05 37.36    15.7,15.84,12.2,12.063,11.475,10.276,
#PN M 1-8   06 53 33.795    +03 08 26.96    11.7,,,,,,

#PN K 3-63
#define the target
target = Target('PN K 3-63', '21 39 11.976', '+55 46 03.94')
#build image stack
stack = Stack(120, 'r-band', 2, 5)
#build image (stack) sequence
sequence = Sequence([stack], 1)
#minimum altitude to observe this target
min_obs_alt = 40.0
#add this observations to the list
observations.append(Observation(observatory, target, sequence, min_obs_alt, user))

#PN K 3-66
#define the target
target = Target('PN K 3-66', '4:36:37.23', '33:39:29.9988')
#build image stack
stack = Stack(120, 'r-band', 2, 5)
#build image (stack) sequence
sequence = Sequence([stack], 1)
#minimum altitude to observe this target
min_obs_alt = 30.0
#add this observations to the list
observations.append(Observation(observatory, target, sequence, min_obs_alt, user))

#IC 2149
#define the target
target = Target('IC 2149', '5:56:23.862', '46:06:17.4996')
#build image stack
stack = Stack(120, 'r-band', 2, 5)
#build image (stack) sequence
sequence = Sequence([stack], 1)
#minimum altitude to observe this target
min_obs_alt = 30.0
#add this observations to the list
observations.append(Observation(observatory, target, sequence, min_obs_alt, user))

#PN M 1-4
#define the target
target = Target('PN M 1-4', '3:41:43.439', '52:16:59.8512')
#build image stack
stack = Stack(120, 'r-band', 2, 5)
#build image (stack) sequence
sequence = Sequence([stack], 1)
#minimum altitude to observe this target
min_obs_alt = 30.0
#add this observations to the list
observations.append(Observation(observatory, target, sequence, min_obs_alt, user))

#PN M 1-6
#define the target
target = Target('PN M 1-6', '06:35:45.126', '-00:05:37.36')
#build image stack
stack = Stack(120, 'r-band', 2, 5)
#build image (stack) sequence
sequence = Sequence([stack], 1)
#minimum altitude to observe this target
min_obs_alt = 30.0
#add this observations to the list
observations.append(Observation(observatory, target, sequence, min_obs_alt, user))

#PN M 1-8
#define the target
target = Target('PN M 1-8', '6:53:33.7949', '3:08:26.9592')
#build image stack
stack = Stack(120, 'r-band', 2, 5)
#build image (stack) sequence
sequence = Sequence([stack], 1)
#minimum altitude to observe this target
min_obs_alt = 30.0
#add this observations to the list
observations.append(Observation(observatory, target, sequence, min_obs_alt, user))

telescope.slackdebug('Starting nebulizer...')

#start up the scheduler
scheduler = Scheduler(observatory, observations)

#loop through the observations
next_observation = scheduler.whatsNext()   

for obs in observations:
    print obs.toString()
sys.exit(1)

count = 0 
while next_observation != None:
       
    #check sun, clouds, slit, etc.
    telescope.checkSun()
    if telescope.is_cracked:
        telescope.checkSlit()
    telescope.checkAlt()
    telescope.checkClouds()       
       
    if Time.now() >= next_observation.min_obs_alt_time: #its time to look at this object
        count += 1
        if not telescope.is_cracked:
            telescope.crackit()
        #perform observation
        telescope.slackdebug("Observations for %s are starting..."%(next_observation.target.getName()))
        #point the scope
        telescope.pinpoint(next_observation)
        #get the images
        telescope.getImage(next_observation)
        telescope.slackdebug("Observations for %s are complete."%(next_observation.target.getName()))
        scheduler.isDone(next_observation) #mark as complete
    else:
        telescope.slackdebug('The next observation (%s) will start in %d min (at %s)...'%(next_observation.target.getName(), (next_observation.min_obs_alt_time-Time.now()).sec/60, next_observation.min_obs_alt_time.iso[:-7]))
        time.sleep(delay_time)
    
    #if its been a while, check our target list again
    next_observation = scheduler.whatsNext()

telescope.slackdebug('Nebulization complete. Performed %d observations.'%count)
