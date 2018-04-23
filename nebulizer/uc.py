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
user='uc'
        
#
#main
#

#set up logger
logger = log.setup_custom_logger('uc')

#simulate? set to True
simulate = False

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

##5358 (1992 QH)
##define the target
#target = Target('5358', '07 38 46.48', '+36 07 12.2')
##build image stack
#stack = Stack(120, 'clear', 2, 1)
##build image (stack) sequence
#sequence = Sequence([stack], Sequence.CONTINUOUS)
##minimum altitude to observe this target
#min_obs_alt = 30.0
##add this observations to the list
#observations.append(Observation(observatory, target, sequence, min_obs_alt, user))
min_obs_alt = 30.0

target = Target('maxgoldberg.M104', '12 39 59.299', '-11 37 21.000')
stack = Stack(240, 'clear', 2, 1)
sequence = Sequence([stack], 1)
observations.append(Observation(observatory, target, sequence, min_obs_alt, user))

target = Target('maxgoldberg.M104', '12 39 59.299', '-11 37 21.000')
stack = Stack(300, 'clear', 2, 1)
sequence = Sequence([stack], 1)
observations.append(Observation(observatory, target, sequence, min_obs_alt, user))

target = Target('shelstrom.M104', '12 39 59.299', '-11 37 21.000')
sequence = Sequence([Stack(60, 'clear', 2, 2), Stack(40, 'r-band', 2, 2), Stack(60, 'g-band', 2, 2)], 1)
observations.append(Observation(observatory, target, sequence, min_obs_alt, user))

target = Target('rich.oddjob.Jupiter', '17 29 06', '-22 37 07')
stack = Stack(0.1, 'z-band', 2, 1)
sequence = Sequence([stack], 1)
observations.append(Observation(observatory, target, sequence, min_obs_alt, user))

target = Target('cepang.pleiades', '03 46 03', '-24 07 57')
stack = Stack(200, 'g-band', 2, 1)
stack2 = Stack(200, 'i-band', 2, 1)
sequence = Sequence([stack, stack2], 1)
observations.append(Observation(observatory, target, sequence, min_obs_alt, user))

target = Target('cepang.05_14_32.-08_12_05', '05 14 32', '-08 12 05')
stack = Stack(200, 'g-band', 2, 1)
stack2 = Stack(200, 'i-band', 2, 1)
sequence = Sequence([stack, stack2], 1)
observations.append(Observation(observatory, target, sequence, min_obs_alt, user))

target = Target('cesmerian.NGC227', '00 42 37', '-01 31 41')
sequence = Sequence([Stack(60, 'clear', 2, 1), Stack(60, 'u-band', 2, 1), Stack(60, 'g-band', 2, 1), Stack(60, 'r-band', 2, 1), Stack(60, 'i-band', 2, 1)], 1)
observations.append(Observation(observatory, target, sequence, min_obs_alt, user))

target = Target('cesmerian.NGC2775', '09 10 20', '07 02 14')
sequence = Sequence([Stack(60, 'clear', 2, 1), Stack(60, 'u-band', 2, 1), Stack(60, 'r-band', 2, 1), Stack(60, 'i-band', 2, 1)], 1)
observations.append(Observation(observatory, target, sequence, min_obs_alt, user))

target = Target('cesmerian.NGC4475', '12 29 48', '27 14 36')
sequence = Sequence([Stack(60, 'clear', 2, 1), Stack(60, 'u-band', 2, 1), Stack(60, 'g-band', 2, 1), Stack(60, 'r-band', 2, 1), Stack(60, 'i-band', 2, 1)], 1)
observations.append(Observation(observatory, target, sequence, min_obs_alt, user))

target = Target('sophiavlahakis.M42', '05 35 17', '-05 23 25')
sequence = Sequence([Stack(60, 'g-band', 2, 1), Stack(60, 'r-band', 2, 1), Stack(60, 'i-band', 2, 1)], 1)
observations.append(Observation(observatory, target, sequence, min_obs_alt, user))

target = Target('sophiavlahakis.M83', '13 37 00', '-29 52 02')
sequence = Sequence([Stack(40, 'clear', 2, 5), Stack(40, 'r-band', 2, 5), Stack(40, 'h-alpha', 2, 5)], 1)
observations.append(Observation(observatory, target, sequence, min_obs_alt, user))

target = Target('mroussi.M42', '05 35 17', '-05 23 25')
sequence = Sequence([Stack(30, 'g-band', 2, 1), Stack(30, 'r-band', 2, 1), Stack(30, 'i-band', 2, 1)], 1)
observations.append(Observation(observatory, target, sequence, min_obs_alt, user))

target = Target('oliviapalid.Neptune', '23 10 31', '-06 21 03')
sequence = Sequence([Stack(300, 'g-band', 2, 2), Stack(300, 'i-band', 2, 2)], 1)
observations.append(Observation(observatory, target, sequence, min_obs_alt, user))

telescope.slackdebug('Starting nebulizer...')

#start up the scheduler
scheduler = Scheduler(observatory, observations)

#loop through the observations
next_observation = scheduler.whatsNext() 

#for obs in observations:
#    print obs.toString()

if next_observation:
    print next_observation.toString()

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
        time.sleep(delay_time)
    
    #if its been a while, check our target list again
    next_observation = scheduler.whatsNext()

telescope.slackdebug('Nebulization complete. Performed %d observations.'%count)
