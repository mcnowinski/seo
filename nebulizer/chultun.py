#
#chultun module
#
#A chultun is a bottle-shaped cavity, excavated by the ancient Maya
#into the softlimestone bedrock typical of the Maya area in the Yucatan peninsula.
#Archaeologists and historians report that chultuns were used for storage purposes,
#for rainwater or other things, and after abandonment for trash and sometimes even burials.
#

import sys
import numpy as np
import astropy.units as u
from astropy.time import Time
from astroplan import Observer
from astropy.coordinates import SkyCoord, EarthLocation, AltAz, get_sun
from astroplan import download_IERS_A
import matplotlib.pyplot as plt
import datetime
import log
import traceback
import log
import subprocess
import time
import re

#set up logger
logger = log.setup_custom_logger('chultun')

#
#the SEO elescope
#
class Telescope():

    #simualte telescope operation?
    simulate = False

    #max. cloud cover, 0-1
    max_clouds_image = 0.4 #cloudier than this will pause imaging
    max_clouds_slit = 0.8 #cloudier than this will close the observatory

    #max sun altitude
    max_sun_alt = -10
    
    #is dome open?
    is_cracked = False

    #min target elevation
    min_alt = 28.0 #no pointing below this elevation    

    def __init__(self, simulate):
        self.simulate = simulate
        
    def done():
        self.squeezeit()
        is_cracked = False
        sys.exit(1)
    
    def runSubprocess(self, command_array, simulate=False, communicate=True):
        #command array is array with command and all required parameters
        if simulate:
            logger.debug('Simulating subprocess "%s".'%(command_array))
            return ('', 0, 0)         
        try:
            sp = subprocess.Popen(command_array, stderr=subprocess.PIPE, stdout=subprocess.PIPE)
            logger.info('Running subprocess ("%s" %s)...'%(' '.join(command_array), sp.pid))
            sp.wait()
            if communicate:
                output, error = sp.communicate(b'\n\n') 
                if error:
                    logger.error(error)            
                return (output, error, sp.pid)
            else: #some processes, like keepopen, hang forever with .communicate()
                return ('', '', 0)
        except Exception as e:
            logger.error(traceback.format_exception(*sys.exc_info()))
            return ('', 'Unknown error.', 0)    
    
    #send alert message to slack        
    def slackalert(self, msg):
        if self.simulate:
            self.slackdev(msg)
            return
        msg = datetime.datetime.utcnow().strftime('%m-%d-%Y %H:%M:%S ') + msg
        logger.debug(msg)
        (output, error, pid) = self.runSubprocess(['slackalert', msg])
            
    #send debug message to slack        
    def slackdebug(self, msg):
        if self.simulate:
            self.slackdev(msg)
            return
        msg = datetime.datetime.utcnow().strftime('%m-%d-%Y %H:%M:%S ') + msg
        logger.debug(msg)
        (output, error, pid) = self.runSubprocess(['slackdebug', msg])

    #send dev message to slack        
    def slackdev(self, msg):
        msg = datetime.datetime.utcnow().strftime('%m-%d-%Y %H:%M:%S ') + msg
        logger.debug(msg)
        (output, error, pid) = self.runSubprocess(['slackdev', msg])    
        
    #send preview of fits image to Slack    
    def slackpreview(self, fits):
        if self.simulate:
            self.slackdev('Placeholder for image (%s).'%fits)
            return
        (output, error, pid) = self.runSubprocess(['stiffy', fits, 'image.tif'])
        (output, error, pid) = self.runSubprocess(['convert', '-resize','50%', '-normalize', '-quality', '75', 'image.tif', 'image.jpg'])
        (output, error, pid) = self.runSubprocess(['slackpreview', 'image.jpg', fits])
        (output, error, pid) = self.runSubprocess(['rm', 'image.jpg', 'image.tif'])

    def squeezeit(self):
        logger.info('Closing the observatory...')
        (output, error, pid) = self.runSubprocess(['squeezeit'], self.simulate) 
        #sys.exit(1)
  
    def crackit(self):
        logger.info('Opening the observatory...')
        (output, error, pid) = self.runSubprocess(['tx', 'lock', 'clear'], self.simulate)
        (output, error, pid) = self.runSubprocess(['tx', 'lock', 'user=mcnowinski', 'email=mcnowinski@gmail.com', 'phone=7032869140', 'comment=chultun'], self.simulate)
        (output, error, pid) = self.runSubprocess(['tin', 'interrupt'], self.simulate)
        (output, error, pid) = self.runSubprocess(['openup', 'nocloud'], self.simulate)
        (output, error, pid) = self.runSubprocess(['keepopen', 'kill', 'debug', 'maxtime=36000', 'slit'], self.simulate, False)
        (output, error, pid) = self.runSubprocess(['tx', 'track', 'on'], self.simulate)        
        (output, error, pid) = self.runSubprocess(['tx', 'slit'])
        match = re.search('slit=open', output)
        if match:
            logger.info('Observatory is open.')
            self.is_cracked = True
        else:
            logger.error('Observatory could not be opened.')        
        #sys.exit(1)
        
    #check slit
    #if the slit is closed, alert the observer via Slack
    def checkSlit(self):
        (output, error, pid) = self.runSubprocess(['tx', 'slit'])
        match = re.search('slit=open', output)
        if match:
            logger.debug('Slit is open.')
        elif self.simulate:
            logger.debug('Slit is (not really) open.')        
        else:
            logger.error('Slit has closed unexpectedly. Trying to re-open it...')
            self.slackalert('Slit has closed unexpectedly. Trying to re-open it...')
            #try to re-open it
            (output, error, pid) = self.runSubprocess(['tx', 'slit', 'open'])
            match = re.search('slit=open', output)
            if match:
                logger.debug('Slit is open.')
            elif self.simulate:
                logger.debug('Slit is (not really) open.')        
            else:
                logger.error('Could not re-open the slit.')
                self.slackalert('Could not re-open the slit.')
                self.done()
            
    #check sun altitude
    #if the sun is too high, squeeze it
    def checkSun(self):
        (output, error, pid) = self.runSubprocess(['sun'])
        match = re.search('alt=([\\-\\+\\.0-9]+)', output)
        if match:
            alt = float(match.group(1))
            logger.debug('Sun altitude is %s deg.'%alt)
            if alt > self.max_sun_alt:
                logger.error('Sun is too high (%s > %s deg).'%(alt, self.max_sun_alt))
                self.slackalert('Sun is too high (%s > %s deg).'%(alt, self.max_sun_alt))
                if not self.simulate:
                    self.done()
        else:
            logger.error('Error. Could not determine the current altitude of the sun (%s).'%output)
            self.slackalert('Error. Could not determine the current altitude of the sun.')
            self.done()
     
    #check altitude of the telescope
    def checkAlt(self):
        (output, error, pid) = self.runSubprocess(['tx','where'])
        match = re.search('alt=([\\-\\+\\.0-9]+)', output)
        if match:
            alt = float(match.group(1))
            logger.debug('Telescope altitude is %s deg.'%alt)
            if alt < self.min_alt:
                logger.error('Telescope altitude is too low (%s < %s deg).'%(alt, self.min_alt))
                self.slackalert('Telescope altitude is too low (%s < %s deg).'%(alt, self.min_alt))
                if not self.simulate:
                    self.done()
        else:
            logger.error('Error. Could not determine the current altitude of the telescope (%s).'%output)
            self.slackalert('Error. Could not determine the current altitude of the telescope.')
            self.done()
        
    #check clouds
    #if its too cloudy, wait it out...
    def checkClouds(self):
        (output, error, pid) = self.runSubprocess(['tx','taux'])
        match = re.search('cloud=([\\-\\.0-9]+)', output)
        if match:
            clouds = float(match.group(1))
            logger.debug('Cloud cover is %d%%.'%int(clouds*100))
            if clouds >= self.max_clouds_slit:
                logger.error('Too many clouds (%d%%). Aborting image sequence...'%int(clouds*100))
                if not self.simulate:
                    self.done()
            while clouds >= self.max_clouds_image:
                self.slackalert('Too many clouds (%d%%). Pausing image sequence...'%int(clouds*100))
                self.checkSun()
                self.checkSlit()
                if self.simulate:
                    break
                time.sleep(30)
                match = re.search('cloud=([\\-\\.0-9]+)', output)
                if match:
                    clouds = float(match.group(1))
                    if clouds >= self.max_clouds_slit:
                        logger.error('Too many clouds (%d%%). Aborting image sequence...'%int(clouds*100))
                        if not self.simulate:
                            self.squeezeit()
                    logger.debug('Cloud cover is %d%%.'%int(clouds*100))
                else:
                    logger.error('Cloud command failed (%s).'%output)            
        return True

    def pinpoint(self, observation):
        name = observation.target.name
        ra = observation.target.getRa().replace(' ', ':')
        dec = observation.target.getDec().replace(' ', ':')
        
        self.slackdebug('Pointing telescope to %s (RA=%s, DEC=%s)...'%(name, ra, dec))
        
        #point the telescope
        start_point = datetime.datetime.utcnow()
        (output, error, pid) = self.runSubprocess(['pinpoint','%s'%ra, '%s'%dec], self.simulate)
        end_point = datetime.datetime.utcnow()
        
        #calculate pointing time in seconds
        dt_point = (end_point-start_point).total_seconds()
        logger.debug('Pinpointing telescope required %d seconds.'%dt_point)
        
        #check the current telescope position
        (output, error, pid) = self.runSubprocess(['tx','where'])    
        
    def getImage(self, observation, do_pinpoint=True):
        for repeat in range(0, observation.sequence.repeat): #how many times does this imaging sequence get repeated?
            for stack in observation.sequence.stacks: #iterate thru the image sets
                for count in range(0, stack.count): #how many times does this stack get repeated?
                    name = observation.target.name
                    filter = stack.filter
                    exposure = stack.exposure
                    binning = stack.binning
                    user = observation.observer
                    #get image
                    fits = '%s_%s_%dsec_bin%d_%s_%s_num%d_seo.fits'%(name, filter, exposure, binning, user, datetime.datetime.utcnow().strftime('%Y%b%d_%Hh%Mm%Ss'), count)
                    fits = fits.replace(' ', '_')
                    fits = fits.replace('(', '')    
                    fits = fits.replace(')', '')
                    self.slackdebug('Taking image (%s)...'%(fits))
                    if self.simulate:
                        #(output, error, pid) = runSubprocess(['image','dark','time=%d'%t_exposure,'bin=%d'%bin, 'outfile=%s'%fits], simulate)
                        time.sleep(exposure)
                        error = 0
                    else:
                        (output, error, pid) = self.runSubprocess(['pfilter','%s'%filter])                          
                        (output, error, pid) = self.runSubprocess(['image','time=%d'%exposure,'bin=%d'%binning, 'outfile=%s'%fits])    
                    if not error:
                        self.slackdebug('Got image (%s).'%fits)
                        if not self.simulate:
                            self.slackpreview(fits)
                        #IMAGE_PATHNAME=$STARS_IMAGE_PATH/`date -u +"%Y"`/`date -u +"%Y-%m-%d"`/${NAME}
                        #(ssh -q -i $STARS_PRIVATE_KEY_PATH $STARS_USERNAME@$STARS_SERVER "mkdir -p $IMAGE_PATHNAME"; scp -q -i $STARS_PRIVATE_KEY_PATH $IMAGE_FILENAME $STARS_USERNAME@$STARS_SERVER:$IMAGE_PATHNAME/$IMAGE_FILENAME) &
                        #(output, error, pid) = runSubprocess(['tostars','%s'%name.replace(' ', '_').replace('(', '').replace(')', ''),'%s'%fits])         
                    else:
                        self.slackdebug('Error. Image command failed (%s).'%fits)
                    if do_pinpoint:
                        self.pinpoint(observation)

#
#an astronomical observation
#
class Observation():

    sequence = None #image sequence
    target = None #target
    observatory = None #observatory
    min_obs_alt = 30 #min alt to start observations in deg
    obs_start_time = None #when will target best be observable?
    min_obs_time = None #when is the target first observable?
    max_obs_time = None #when is the target last observable?
    active = True #is this observation (still) active?
    id = -1 #id
    observer = None
    
    #init
    def __init__(self, observatory, target, sequence, min_obs_alt, observer):
        self.observatory = observatory
        self.target = target
        self.sequence = sequence
        self.min_obs_alt = min_obs_alt
        self.observer = observer
    
    def toString(self):
        return '%s\n%s\n%smin_alt=%f deg\nobs_time=%s\nid=%d\nactive=%d\nmin_obs_time=%s\nmax_obs_time=%s'%(self.observatory.toString(), self.target.toString(), self.sequence.toString(), self.min_obs_alt, self.obs_start_time, self.id, self.active, self.min_obs_time, self.max_obs_time)
        
    def setObservationStartTime(self, start_time):
        self.obs_start_time = start_time

    def setObservationMinTime(self, min_time):
        self.min_obs_time = min_time        

    def setObservationMaxTime(self, max_time):
        self.max_obs_time = max_time  
        
#
#the astronomical target
#
class Target():
   
    #init
    def __init__(self, name, ra, dec):
        self.name = name
        self.ra = ra #hour:min:sec
        self.dec = dec #deg:min:sec

    #name   
    def getName(self):
        return self.name
        
    def setName(self, name):
        self.name = name        
        
    #ra = right ascension
    #eventually expand to allow current coord lookup based on name?    
    def getRa(self):
        return self.ra
        
    def setRa(self, ra):
        self.ra = ra

    #dec = declination
    #eventually expand to allow current coord lookup based on name?      
    def getDec(self):
        return self.dec
        
    def setDec(self, dec):
        self.dec = dec
        
    def toString(self):
        return 'target: name=%s, ra=%s, dec=%s'%(self.name, self.ra, self.dec)

#
#settings for a single set of astronomical images
#
class Stack():

    exposure = 10 #exposure time in seconds
    filter = 'clear' #filter, e.g., clear, h-alpha, u-band, g-band, r-band, i-band, z-band
    binning = 1 #binning, e.g. 1 or 2
    count = 1 #number of images in this stack
   
    #init
    def __init__(self, exposure, filter, binning, count):
        self.exposure = exposure
        self.filter = filter
        self.binning = binning
        self.count = count
    
    def toString(self):
        return 'image stack: exposure=%d, filter=%s, binning=%d, count=%d'%(self.exposure, self.filter, self.binning, self.count)

#
#sequence of astronomical image stacks
#
class Sequence():

    stacks = [] #list of image stacks
    repeat = 1 #number of times to repeat this sequence
    
    #init
    def __init__(self, stacks, repeat):
        self.stacks = stacks
        self.repeat = repeat
    
    def addStack(self, stack):
        self.stacks.append(stack)
    
    def toString(self):
        sequence_string = 'sequence: repeat=%d\n'%(self.repeat)
        for stack in self.stacks:
            sequence_string += '  %s\n'%stack.toString()
        return sequence_string

#
#the observatory
#        
class Observatory():

    code = None
    latitude = 0.0 #in decimal degrees
    longitude = 0.0 #in decimal degrees
    altitude = 0.0 #in meters
    timzeone = None
    
    #init
    def __init__(self, code, latitude, longitude, altitude, timezone):
        self.code = code
        self.latitude = latitude    
        self.longitude = longitude  
        self.altitude = altitude 
        self.timezone = timezone

    def toString(self):
        observatory_string = 'observatory: code=%s, lat=%f, lon=%f, alt=%f'%(self.code, self.latitude, self.longitude, self.altitude)
        return observatory_string
        

class Scheduler():
#based on Amanda Pagul's schedule function for the SEO queue
# This function receives a list of target coordinates and outputs a primary
# observable target
# Inputs:
# -------
# target_list :str: :list: List containing the names of the targets (i.e. [(id, ra, dec), (id, ra, dec), ...]).
# endtime: :list: 'obj' datetime (i.e. datetime(year, month, day, hour, minute, second))
# Outputs:
# --------
# primary_target :str: (id, ra, dec) of the highest priority target
# wait :astropy.Time: Wait time in seconds until optimal observation.
# It takes the value -1 when the object(s) is not observable.
#
    observatory = None #an Observatory
    observations = None #list of Observations
    last_id = 0 #unique id for each observation
    
    sunset_time = None #nearest sunset
    sunrise_time = None #next sunrise
    
    #max_sun_alt = -12 #what defines dark? (deg)
    #min_target_alt = 30 #how low can you go? (deg)
    
    def __init__(self, observatory, observations):
        self.observatory = observatory
        self.observations = observations
        #assign unique ids
        for observation in self.observations:
            observation.id = self.last_id
            self.last_id += 1
        
        #uncomment to download the latest for astroplan
        logger.debug('Updating Astroplan IERS Bulletin A...')
        download_IERS_A()
        
        #get *nearest* sunset and *next* sunrise times
        #still not a big fan of this!
        observatory_location_obsplan = Observer(longitude=self.observatory.longitude*u.deg, latitude=self.observatory.latitude*u.deg, elevation=self.observatory.altitude*u.m, name=self.observatory.code, timezone=self.observatory.timezone) 
        self.sunset_time = observatory_location_obsplan.twilight_evening_nautical(Time.now(), which="nearest") 
        self.sunrise_time = observatory_location_obsplan.twilight_morning_nautical(Time.now(), which="next") 
        logger.debug('The nearest sunset is %s. The next sunrise is %s.'%(self.sunset_time.iso, self.sunrise_time.iso))       

    def isDone(self, observation):
        for obs in self.observations:
            if obs.id == observation.id:
                obs.active = False
                return
        logger.error('Matching observation (%d) not found.'%observation.id)
        
    def whatsNext(self):   
        obs_stats = {'target':[], 'altitude':[], 'time':[], 'wait':[], 'starttime':[], 'id':[]}   
        
        observatory_location = EarthLocation(lat=self.observatory.latitude*u.deg, lon=self.observatory.longitude*u.deg, height=self.observatory.altitude*u.m)            
        
        #build alt-az coordinate frame for observatory over next 15 hours
        #obs_time = Time.now()
        obs_time = self.sunset_time
        delta_obs_time = np.linspace(0, (self.sunrise_time-obs_time).sec/3600., 1000)*u.hour
        #the next 15 hours
        times = obs_time+delta_obs_time
        #build frame
        frame = AltAz(obstime=times, location=observatory_location)
        #print frame
                
        #loop thru observations, suggest the next best target
        for observation in self.observations:
                
            if observation.active == False:
                logger.info('This observation (%s) is no (longer) active. Skipping...'%observation.target.getName())
                continue
                
            target_ra = observation.target.getRa()
            target_dec = observation.target.getDec()
            input_coordinates = target_ra+" "+target_dec
            
            obs_stats['target'].append(observation.target.getName())

            try:
                target_coordinates = SkyCoord(input_coordinates, unit=(u.hourangle, u.deg))
            except:
                continue
            
            target_altaz = target_coordinates.transform_to(frame)
            
            #if (np.max(target_altaz.alt)) > observation.min_obs_alt*u.degree:
            #    #obs_stats['altitude'].append(np.max(target_altaz.alt))
            #    #get time of min observable altitude
            #    min_alt_time = times[np.where(target_altaz.alt > observation.min_obs_alt*u.degree)] 
            #    observation.setObservationStartTime(min_alt_time[0])
            ##else:
            ##    obs_stats['altitude'].append(0*u.degree)
                       
            max_alt_time = times[np.argmax(target_altaz.alt)]
            obs_stats['time'].append(max_alt_time)
                       
            #print max_alt_time.iso
            
            aux_delta_time = delta_obs_time[np.argmax(target_altaz.alt)]
        
            #print obs_stats['altitude'][i], times[np.argmax(target_altaz.alt)], sunset_time, sunrise_time
            if ((times[np.argmax(target_altaz.alt)] >= self.sunset_time) & (times[np.argmax(target_altaz.alt)] < self.sunrise_time)):
            #if (obs_stats['altitude'][observation.id]>0*u.degree) & (times[np.argmax(target_altaz.alt)] > sunset_time) & (times[np.argmax(target_altaz.alt)] < sunrise_time):
                obs_stats['wait'].append(aux_delta_time.to(u.second))
            else:
                obs_stats['wait'].append(-1*u.s)
             
            #track id of this observation
            obs_stats['id'].append(observation.id)
            
            #get min and max obs times
            #obs_stats['min_obs_time'].append(np.min(times[np.where(target_altaz.alt > observation.min_obs_alt*u.degree)]))
            #obs_stats['max_obs_time'].append(np.max(times[np.where(target_altaz.alt > observation.min_obs_alt*u.degree)]))
            ##set min and max obs times
            self.observations[observation.id].setObservationMinTime(Time(np.min(times[np.where(target_altaz.alt > observation.min_obs_alt*u.degree)])))
            self.observations[observation.id].setObservationMaxTime(Time(np.max(times[np.where(target_altaz.alt > observation.min_obs_alt*u.degree)])))

                
        #print obs_stats['time']
        #print obs_stats['target']        
        #print obs_stats['wait'] 
                
        obs_stats['time']=np.array(obs_stats['time'])
        #print obs_stats['wait']
        #sys.exit(1)
        good_object = np.array([obs_stats['wait'][itgt]>-1*u.s for itgt in range(len(obs_stats['wait']))]) 
        
        if np.count_nonzero(good_object)>0:
            if np.count_nonzero(good_object)>1:
                aux_id = np.argmin(Time(obs_stats['time'][good_object])-Time.now())
                primary_target_id = np.where(good_object)[0][aux_id]
                #print obs_stats['target'][primary_target_id]
                primary_target = np.array(obs_stats['target'])[primary_target_id]
                #print primary_target_id
            else:
                primary_target_id = np.where(good_object)[0][0] 
                primary_target = np.array(obs_stats['target'])[primary_target_id]
            #set start time
            self.observations[obs_stats['id'][primary_target_id]].setObservationStartTime(Time(obs_stats['time'][primary_target_id]))
            #logger.debug('target_id=%d'%obs_stats['id'][primary_target_id])            
        else:
            return None
                  
        return self.observations[obs_stats['id'][primary_target_id]]
