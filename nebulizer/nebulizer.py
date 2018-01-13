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
simulate = True

#time between checks for object observability in seconds
delay_time = 30

#min time available for background observations in seconds
min_background_time = 10*60 #10 minutes

#

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

# #1639 Rollandia
# #define the target
# target = Target('1639 Rollandia', '04 22 19.65', '+18 20 26.8')
# #build image stack
# stack = Stack(120, 'clear', 2, 1)
# #build image (stack) sequence
# sequence = Sequence([stack], Sequence.CONTINUOUS)
# #minimum altitude to observe this target
# min_obs_alt = 30.0
# #add this observations to the list
# observations.append(Observation(observatory, target, sequence, min_obs_alt, user))

#5358 (1992 QH)
#define the target
target = Target('5358', '07 38 46.48', '+36 07 12.2')
#build image stack
stack = Stack(120, 'clear', 2, 1)
#build image (stack) sequence
sequence = Sequence([stack], Sequence.CONTINUOUS)
#minimum altitude to observe this target
min_obs_alt = 30.0
#add this observations to the list
observations.append(Observation(observatory, target, sequence, min_obs_alt, user))

#PN K 3-63
#define the target
target = Target('PN K 3-63', '21 39 11.976', '+55 46 03.94')
#build image stack
stack = Stack(120, 'r-band', 2, 5)
#build image (stack) sequence
sequence = Sequence([stack], 1)
#minimum altitude to observe this target
min_obs_alt = 30.0
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

in_between_observation = scheduler.whatsInBetween()  

#for obs in observations:
#    print obs.toString()

if next_observation:
    print next_observation.toString()

if in_between_observation:
    print in_between_observation.toString()

#wait for sun to set
telescope.checkSun(True)

count = 0 
while next_observation != None:
       
    #check sun, clouds, slit, etc.
    telescope.checkSun()
    if telescope.is_cracked:
        telescope.checkSlit()
    telescope.checkAlt()
    telescope.checkClouds()       
       
    if Time.now() >= next_observation.min_obs_time: #its time to look at this object
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
        wait_time_s = (next_observation.min_obs_time-Time.now()).sec
        telescope.slackdebug('The next observation (%s) will start in %d min (at %s)...'%(next_observation.target.getName(), wait_time_s/60, next_observation.min_obs_time.iso[:-7]))
        #should we do a background observation while we are waiting?
        if wait_time_s >= min_background_time:
            in_between_observation = scheduler.whatsInBetween()
            if in_between_observation != None:
                if not telescope.is_cracked:
                    telescope.crackit()
                #perform observation
                telescope.slackdebug("Background observations for %s are starting..."%(in_between_observation.target.getName()))
                #point the scope
                telescope.pinpoint(in_between_observation)
                #get the images
                telescope.getImage(in_between_observation, next_observation.min_obs_time)

        time.sleep(delay_time)
    
    #if its been a while, check our target list again
    next_observation = scheduler.whatsNext()

#all foreground observations are complete, work the background observations
in_between_observation = scheduler.whatsInBetween()
while in_between_observation != None:
       
    #check sun, clouds, slit, etc.
    telescope.checkSun()
    if telescope.is_cracked:
        telescope.checkSlit()
    telescope.checkAlt()
    telescope.checkClouds()

    if not telescope.is_cracked:
        telescope.crackit()
    #perform observation
    telescope.slackdebug("Background observations for %s are starting..."%(in_between_observation.target.getName()))
    #point the scope
    telescope.pinpoint(in_between_observation)
    #get the images
    telescope.getImage(in_between_observation, in_between_observation.max_obs_time)

    in_between_observation = scheduler.whatsInBetween()

telescope.slackdebug('Nebulization complete. Performed %d observations.'%count)
