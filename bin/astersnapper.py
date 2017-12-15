import callhorizons
import datetime
import os
import sys
import subprocess
import time
import logging
import unicodedata
import string
import re
import math

from astropy import units as u
from astropy.coordinates import Angle
from logging.handlers import RotatingFileHandler
from dateutil import parser
from astropy.time import Time

#run external process; track output, errors, and pid
def runSubprocess(command_array, simulate=False):
    #command array is array with command and all required parameters
    if simulate:
        logger.debug('Simulating subprocess "%s".'%(command_array))
        return ('', 0, 0)    
    try:
        sp = subprocess.Popen(command_array, stderr=subprocess.PIPE, stdout=subprocess.PIPE)
        logger.info('Running subprocess "%s" (%s)...'%(' '.join(command_array), sp.pid))
        output, error = sp.communicate()
        logger.debug(output)
        if error:
            logger.error(error)
        return (output, error, sp.pid)
    except:
        logger.error('Error. Subprocess ("%s") failed.'%(' '.join(command_array)))
        return ('', 'Unknown error.', 0)     

#send alert message to slack        
def slackalert(msg):
    if simulate:
        slackdev(msg)
        return
    msg = datetime.datetime.utcnow().strftime('%m-%d-%Y %H:%M:%S ') + msg
    logger.debug(msg)
    (output, error, pid) = runSubprocess(['slackalert', msg])
        
#send debug message to slack        
def slackdebug(msg):
    if simulate:
        slackdev(msg)
        return
    msg = datetime.datetime.utcnow().strftime('%m-%d-%Y %H:%M:%S ') + msg
    logger.debug(msg)
    (output, error, pid) = runSubprocess(['slackdebug', msg])

#send dev message to slack        
def slackdev(msg):
    msg = datetime.datetime.utcnow().strftime('%m-%d-%Y %H:%M:%S ') + msg
    logger.debug(msg)
    (output, error, pid) = runSubprocess(['slackdev', msg])    
    
#send preview of fits image to Slack    
def slackpreview(fits):
    if simulate:
        slackdev('Placeholder for image (%s).'%fits)
        return
    (output, error, pid) = runSubprocess(['stiffy', fits, 'image.tif'])
    (output, error, pid) = runSubprocess(['convert', '-resize','50%', '-normalize', '-quality', '75', 'image.tif', 'image.jpg'])
    (output, error, pid) = runSubprocess(['slackpreview', 'image.jpg', fits])
    (output, error, pid) = runSubprocess(['rm', 'image.jpg', 'image.tif'])

def squeezeit():
    logger.info('Closing the observatory...')
    (output, error, pid) = runSubprocess(['squeezeit'], simulate) 
    sys.exit(1)
    
#check slit
#if the slit is closed, alert the observer via Slack
def checkSlit():
    (output, error, pid) = runSubprocess(['tx', 'slit'])
    match = re.search('slit=open', output)
    if match:
        logger.debug('Slit is open.')
    elif simulate:
        logger.debug('Slit is (not really) open.')        
    else:
        #send a repeated alert to Slack
        while True:
            slackalert('Slit has closed unexpectedly.')
            time.sleep(20)
    return True

#check sun altitude
#if the sun is too high, squeeze it
def checkSun():
    (output, error, pid) = runSubprocess(['sun'])
    match = re.search('alt=([\\-\\+\\.0-9]+)', output)
    if match:
        alt = float(match.group(1))
        logger.debug('Sun altitude is %s deg.'%alt)
        if alt > max_sun_alt:
            logger.info('Sun is too high (%s > %s deg).'%(alt, max_sun_alt))
            if not simulate:
                squeezeit()
    else:
        logger.error('Error. Could not determine the current altitude of the sun (%s).'%output)
   
#check clouds
#if its too cloudy, wait it out...
def checkClouds():
    (output, error, pid) = runSubprocess(['tx','taux'])
    match = re.search('cloud=([\\-\\.0-9]+)', output)
    if match:
        clouds = float(match.group(1))
        logger.debug('Cloud cover is %d%%.'%int(clouds*100))
        if clouds >= max_clouds_slit:
            logger.error('Too many clouds (%d%%). Aborting image sequence...'%int(clouds*100))
            if not simulate:
                squeezeit()
        while clouds >= max_clouds_image:
            slackalert('Too many clouds (%d%%). Pausing image sequence...'%int(clouds*100))
            checkSun()
            if simulate:
                break
            time.sleep(30)
            match = re.search('cloud=([\\-\\.0-9]+)', output)
            if match:
                clouds = float(match.group(1))
                if clouds >= max_clouds_slit:
                    logger.error('Too many clouds (%d%%). Aborting image sequence...'%int(clouds*100))
                    if not simulate:
                        squeezeit()
                logger.debug('Cloud cover is %d%%.'%int(clouds*100))
            else:
                logger.error('Cloud command failed (%s).'%output)            
    return True

#
# main
#    
 
# 
#CHANGE AS NEEDED!
#Some should eventually move to command line parameters or a config file... 
#
#the target
target = '1637'

#move center of FOV (e.g., to avoid bright star)
dRA = 0.0
dDEC = 0.0

#the observatory
observatory_code = 'G52'
observatory_lat = 38.259 #deg
observatory_lon = -122.440 #deg
observatory_elev = 63.8 #m
observatory_tz = 'US/Pacific' #for pytz

#ccd min dimension in arcsec
ccd_min_dim = 25.6*60 #25.6 arcmin
ccd_min_dim_margin = 0.1 #0-1, how much of the ccd_min_dim will we use

#exposure time in seconds
t_exposure = 15

#filter
filter = 'clear'

#binning
bin=2

#time alloted to perform telescope (pin)pointing in seconds
t_pointing = 40

#user, hardcode for now
user='mcnowinski'

#max. cloud cover, 0-1
max_clouds_image = 0.4 #cloudier than this will pause imaging
max_clouds_slit = 0.8 #cloudier than this will close the observatory

#min target elevation
min_alt = 28.0 #no pointing below this elevation

#simulate? set to True
simulate = True

#max sun altitude
max_sun_alt=-10 #no imaging above this sun altitude

#callhorizons constants
#max airmass
max_airmass = 2.0 #30 deg elevation

#configure logging
log_file='/home/mcnowinski/var/log/astersnapper.log'
logger = logging.getLogger('astersnapper')
logging.Formatter.converter = time.gmtime
logger.setLevel(logging.DEBUG)
handler = RotatingFileHandler(log_file, maxBytes=5*1024*1024, backupCount=10)
#handler.setFormatter(logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s'))
handler.setFormatter(logging.Formatter('%(asctime)s %(levelname)s\t%(message)s'))
logger.addHandler(handler)
logger.info('Starting astersnapper...')

if simulate:
    logger.info('Running in simulate mode...')    

#where is the asteroid going to be at the *midpoint* of the exposure?
#first, calculate the *time* of the *midpoint* of the exposure: t = now + t_pointing (in sec) + t_exposure/2 (in sec)
#second, calculate the RA/DEC of the target at *that* time using JPL Horizons (via callhorizons)
start = datetime.datetime.utcnow()
#start += datetime.timedelta(seconds=t_pointing)
#start += datetime.timedelta(seconds=t_exposure/2.0)

#JPL Horizons requires start *and* end times where end > start (by at least 1 minute!)
#end = start + datetime.timedelta(seconds=60)
end = start + datetime.timedelta(days=1)

#get ephemerides for target in JPL Horizons from start to end times
# +------------------+-----------------------------------------------+
# | Property         | Definition                                    |
# +==================+===============================================+
# | targetname       | official number, name, designation [string]   |
# +------------------+-----------------------------------------------+
# | H                | absolute magnitude in V band (float, mag)     |
# +------------------+-----------------------------------------------+
# | G                | photometric slope parameter (float)           |
# +------------------+-----------------------------------------------+
# | datetime         | epoch date and time (str, YYYY-MM-DD HH:MM:SS)|
# +------------------+-----------------------------------------------+
# | datetime_jd      | epoch Julian Date (float)                     |
# +------------------+-----------------------------------------------+
# | solar_presence   | information on Sun's presence (str)           |
# +------------------+-----------------------------------------------+
# | lunar_presence   | information on Moon's presence (str)          |
# +------------------+-----------------------------------------------+
# | RA               | target RA (float, J2000.0)                    |
# +------------------+-----------------------------------------------+
# | DEC              | target DEC (float, J2000.0)                   |
# +------------------+-----------------------------------------------+
# | RA_rate          | target rate RA*cos(DEC) (float, arcsec/s)     |
# +------------------+-----------------------------------------------+
# | DEC_rate         | target rate DEC (float, arcsec/s)             |
# +------------------+-----------------------------------------------+
# | AZ               | Azimuth meas East(90) of North(0) (float, deg)|
# +------------------+-----------------------------------------------+
# | EL               | Elevation (float, deg)                        |
# +------------------+-----------------------------------------------+
# | airmass          | target optical airmass (float)                |
# +------------------+-----------------------------------------------+
# | magextinct       | V-mag extinction due airmass (float, mag)     |
# +------------------+-----------------------------------------------+
# | V                | V magnitude (comets: total mag) (float, mag)  |
# +------------------+-----------------------------------------------+
# | illumination     | fraction of illuminated disk (float)          |
# +------------------+-----------------------------------------------+
# | EclLon           | heliocentr. ecl. long. (float, deg, J2000.0)  |
# +------------------+-----------------------------------------------+
# | EclLat           | heliocentr. ecl. lat. (float, deg, J2000.0)   |
# +------------------+-----------------------------------------------+
# | ObsEclLon        | obscentr. ecl. long. (float, deg, J2000.0)    |
# +------------------+-----------------------------------------------+
# | ObsEclLat        | obscentr. ecl. lat. (float, deg, J2000.0)     |
# +------------------+-----------------------------------------------+
# | r                | heliocentric distance (float, au)             |
# +------------------+-----------------------------------------------+
# | r_rate           | heliocentric radial rate  (float, km/s)       |
# +------------------+-----------------------------------------------+
# | delta            | distance from the observer (float, au)        |
# +------------------+-----------------------------------------------+
# | delta_rate       | obs-centric radial rate (float, km/s)         |
# +------------------+-----------------------------------------------+
# | lighttime        | one-way light time (float, s)                 |
# +------------------+-----------------------------------------------+
# | elong            | solar elongation (float, deg)                 |
# +------------------+-----------------------------------------------+
# | elongFlag        | app. position relative to Sun (str)           |
# +------------------+-----------------------------------------------+
# | alpha            | solar phase angle (float, deg)                |
# +------------------+-----------------------------------------------+
# | sunTargetPA      | PA of Sun->target vector (float, deg, EoN)    |
# +------------------+-----------------------------------------------+
# | velocityPA       | PA of velocity vector (float, deg, EoN)       |
# +------------------+-----------------------------------------------+
# | GlxLon           | galactic longitude (float, deg)               |
# +------------------+-----------------------------------------------+
# | GlxLat           | galactic latitude  (float, deg)               |
# +------------------+-----------------------------------------------+
# | RA_3sigma        | 3sigma pos. unc. in RA (float, arcsec)        |
# +------------------+-----------------------------------------------+
# | DEC_3sigma       | 3sigma pos. unc. in DEC (float, arcsec)       |
# +------------------+-----------------------------------------------+
ch=callhorizons.query(target.upper(), smallbody=True)
ch.set_epochrange(start.strftime("%Y/%m/%d %H:%M:%S"), end.strftime("%Y/%m/%d %H:%M:%S"), '1m')
ch.get_ephemerides(observatory_code, skip_daylight=True, airmass_lessthan=max_airmass)

#check callhorizons results
mins = len(ch)
if mins > 0: #is there something to look at?
    logger.debug('name=%s, dt=%s, RA=%s, DEC=%s, EL=%s, AZ=%s, RArate=%s, DECrate=%s'%(ch['targetname'][0], ch['datetime'][0], ch['RA'][0], ch['DEC'][0], ch['EL'][0], ch['AZ'][0], ch['RA_rate'][0], ch['DEC_rate'][0]))
else:
    logger.error('Error. Could not obtain ephemerides for target (%s).'%target)
    os.sys.exit(1)
 
#object name
name=ch['targetname'][0]
#speed of asteroid in arcsec/sec    
rate = math.sqrt(float(ch['RA_rate'][0])*float(ch['RA_rate'][0]) + float(ch['DEC_rate'][0])*float(ch['DEC_rate'][0]))
#total distance in arcsec spanned by asteroid over mins
span = rate * mins * 60    
#calculate lst
lst = Time.now().sidereal_time('mean', longitude='%fd'%observatory_lon)
#print lst
#print Angle(float(ch['RA'][0])*u.deg, u.hourangle)
ha = lst-Angle(float(ch['RA'][0])*u.deg, u.hourangle)
#print ha
logger.debug('rate=%s, duration=%d min, span=%f arcsec,  lst=%s,ha=%s'%(rate, mins, span, lst, ha))
 
#get start and end times for observation from callhorizons results
start = parser.parse(ch['datetime'][0])
end = parser.parse(ch['datetime'][mins-1])
logger.debug('start=%s, end=%s'%(start.strftime("%Y-%m-%dT%H:%M:%S"), end.strftime("%Y-%m-%dT%H:%M:%S")))
 
session_stop_mins = mins #when does the current session end?
if (span <= ccd_min_dim*ccd_min_dim_margin): #asteroid will remain in the FOV throughout one observing session
    logger.debug('num_sessions=1, span=%f arcsec, fov=%f arsec, margin=%f'%(span, ccd_min_dim, ccd_min_dim_margin))
    #use midpoint of session for pointing
    ra=Angle(float(ch['RA'][mins/2])*u.deg).to_string(unit=u.hour, sep=':')
    dec=Angle(float(ch['DEC'][mins/2])*u.deg).to_string(unit=u.degree, sep=':')
    #logger.debug('Pointing telescope to RA=%s, DEC=%s.'%(ra, dec))
else: #we will need to perform multiple sessions, with re-pointing in between
    #use midpoint of time it takes for object to traverse FOV
    session_stop_mins = int((ccd_min_dim*ccd_min_dim_margin / rate) / 60)
    session_mid_mins = session_stop_mins/2
    ra=Angle(float(ch['RA'][session_mid_mins])*u.deg).to_string(unit=u.hour, sep=':')
    dec=Angle(float(ch['DEC'][session_mid_mins])*u.deg).to_string(unit=u.degree, sep=':')
    logger.debug('num_sessions=%d, span=%f arcsec, fov=%f arsec, margin=%f'%(mins/session_stop_mins, span, ccd_min_dim, ccd_min_dim_margin))
 
#wait for asteroid to be observable
while ((start-datetime.datetime.utcnow()).total_seconds() > 0):
    logger.info('Waiting for %s to be observable at %s...'%(name, start.strftime("%Y-%m-%dT%H:%M:%S")))
    time.sleep(10)
   
#main loop
count = 0 #image counter
min = 0 #current minute relative to start
sess = 0 #session number
logger.info('Starting observing session %d for %s...'%(sess, name))
while ((end-datetime.datetime.utcnow()).total_seconds() > 0):
    count += 1
    
    #calc current minute relative to start
    min = (datetime.datetime.utcnow()-start).total_seconds()/60

    #logger.debug('Current minute in observing session %d is %d (of %d).'%(min, mins))
    
    if min > session_stop_mins: #object has moved outside FOV, re-point and start new session
        sess += 1
        session_stop_mins = min + int(ccd_min_dim*ccd_min_dim_margin / rate * 60)
        session_mid_mins = min + int(ccd_min_dim*ccd_min_dim_margin / rate * 60)/2
        ra=Angle(float(ch['RA'][session_mid_mins])*u.deg).to_string(unit=u.hour, sep=':')
        dec=Angle(float(ch['DEC'][session_mid_mins])*u.deg).to_string(unit=u.degree, sep=':')   
        logger.info('Starting observing session %d for %s...'%(sess, name))
    
    time.sleep(10)
    
    #ensure tracking is on, tx track on
    (output, error, pid) = runSubprocess(['tx','track','on'], simulate)    
    
    checkSun()
    checkSlit()
    checkClouds()
    ##checkTargetAlt()    
    
    slackdebug('Pointing telescope to %s (RA=%s, DEC=%s, AZ=%s, ALT=%s)...'%(name, ra, dec, ch['AZ'][0], ch['EL'][0]))
    
    #point the telescope
    start_point = datetime.datetime.utcnow()
    (output, error, pid) = runSubprocess(['pinpoint','%s'%ra, '%s'%dec], simulate)
    end_point = datetime.datetime.utcnow()
    
    #calculate pointing time in seconds
    dt_point = (end_point-start_point).total_seconds()
    logger.debug('Pinpointing telescope required %d seconds.'%dt_point)
    
    #check the current telescope position
    (output, error, pid) = runSubprocess(['tx','where'])
    
    #get image
    fits = '%s_%s_%dsec_bin%d_%s_%s_num%d_seo.fits'%(name, filter, t_exposure, bin, user, datetime.datetime.utcnow().strftime('%Y%b%d_%Hh%Mm%Ss'), count)
    fits = fits.replace(' ', '_')
    fits = fits.replace('(', '')    
    fits = fits.replace(')', '')
    slackdebug('Taking image (%s)...'%(fits))
    if simulate:
        (output, error, pid) = runSubprocess(['image','dark','time=%d'%t_exposure,'bin=%d'%bin, 'outfile=%s'%fits], simulate)
        time.sleep(t_exposure)
    else:
        (output, error, pid) = runSubprocess(['image','time=%d'%t_exposure,'bin=%d'%bin, 'outfile=%s'%fits])    
    if not error:
        slackdebug('Got image (%s).'%fits)
        slackpreview(fits)
        #IMAGE_PATHNAME=$STARS_IMAGE_PATH/`date -u +"%Y"`/`date -u +"%Y-%m-%d"`/${NAME}
        #(ssh -q -i $STARS_PRIVATE_KEY_PATH $STARS_USERNAME@$STARS_SERVER "mkdir -p $IMAGE_PATHNAME"; scp -q -i $STARS_PRIVATE_KEY_PATH $IMAGE_FILENAME $STARS_USERNAME@$STARS_SERVER:$IMAGE_PATHNAME/$IMAGE_FILENAME) &
        #(output, error, pid) = runSubprocess(['tostars','%s'%name.replace(' ', '_').replace('(', '').replace(')', ''),'%s'%fits])         
    else:
        slackdebug('Error. Image command failed (%s).'%fits) 
    ##time.sleep(t_exposure+5)
    
    ##calc new position
    #start = datetime.datetime.utcnow()
    #start += datetime.timedelta(seconds=t_pointing)
    #start += datetime.timedelta(seconds=t_exposure/2.0)
    #end = start + datetime.timedelta(seconds=60)
    #slackdebug('Calculating new position for %s (%s)...'%(name, start.strftime("%Y/%m/%d %H:%M:%S")))
    #ch.set_epochrange(start.strftime("%Y/%m/%d %H:%M:%S"), end.strftime("%Y/%m/%d %H:%M:%S"), '1m')
    #ch.get_ephemerides(observatory)
    #logger.debug('name=%s,dt=%s,RA=%s,DEC=%s,EL=%s,AZ=%s'%(name, ch['datetime'][0], ch['RA'][0], ch['DEC'][0], ch['EL'][0], ch['AZ'][0]))
    
logger.info('Stopping astersnapper...')

#close up shop
(output, error, pid) = runSubprocess(['squeezeit'], simulate)     
    