#
# chultun module
#
# A chultun is a bottle-shaped cavity, excavated by the ancient Maya
# into the softlimestone bedrock typical of the Maya area in the Yucatan peninsula.
# Archaeologists and historians report that chultuns were used for storage purposes,
# for rainwater or other things, and after abandonment for trash and sometimes even burials.
#

import sys
import os
import numpy as np
import astropy.units as u
from astropy.time import Time
from astroplan import Observer
from astropy.coordinates import SkyCoord, EarthLocation, AltAz, get_sun, Angle
from astroplan import download_IERS_A
from astropy.io.fits import getheader
import datetime
import log
import traceback
import log
import subprocess
import time
import re
from astroquery.simbad import Simbad
import warnings
import pathlib2
import urllib2
import ch
import glob2
import json

# set up logger
logger = log.get_logger('chultun.log')

# set up dumper
dumper = log.get_dumper('chultun.dmp')

# python path
python_path = '/home/mcnowinski/anaconda2/bin/python'

# path to pinpoint script
pinpoint_path = '/home/mcnowinski/seo/bin/pinpoint.py'

# image path
image_path = '/home/mcnowinski/itzamna/images'

# image archive path
image_archive_path = '/home/mcnowinski/itzamna/archive'

# path to solve-field astrometry executable
solve_field_path = '/home/mcnowinski/astrometry/bin/solve-field'

# path to secret json config file
secret_json_field_path = '/home/mcnowinski/seo/nebulizer/secret.json'

# path to ssh
ssh_path = 'ssh'

# path to scp
scp_path = 'scp'

# max time to allow for (just) pinpoint operation
# report alert if exceeded
max_pinpoint_time_s = 60

#
# the SEO telescope
#


class Telescope():

    # simualte telescope operation?
    simulate = False

    # max. cloud cover, 0-1
    max_clouds_image = 0.7  # cloudier than this will pause imaging
    max_clouds_slit = 0.9  # cloudier than this will close the observatory

    # max sun altitude
    max_sun_alt = -10

    # is dome open?
    is_cracked = False

    # min target elevation
    min_alt = 25.0  # no pointing below this elevation

    # secrets
    secrets = []

    def __init__(self, simulate):
        self.simulate = simulate
        with open(secret_json_field_path, "r") as f:
            self.secrets = json.load(f)

    def done(self):
        self.squeezeit()
        is_cracked = False
        sys.exit(0)

    def runSubprocess(self, command_array, simulate=False, communicate=True):
        # command array is array with command and all required parameters
        if simulate:
            logger.debug('Simulating subprocess "%s".' % (command_array))
            return ('', 0, 0)
        try:
            sp = subprocess.Popen(
                command_array, stderr=subprocess.PIPE, stdout=subprocess.PIPE)
            logger.info('Running subprocess ("%s" %s)...' %
                        (' '.join(command_array), sp.pid))
            sp.wait()
            if communicate:
                output, error = sp.communicate(b'\n\n')
                if error:
                    logger.error('Process (%s) reported error.' %
                                 (command_array))
                    dumper.error(error)
                return (output, error, sp.pid)
            else:  # some processes, like keepopen, hang forever with .communicate()
                return ('', '', 0)
        except Exception as e:
            logger.error(traceback.format_exception(*sys.exc_info()))
            return ('', 'Unknown error.', 0)

    # send alert message to slack
    def slackalert(self, msg):
        if self.simulate:
            self.slackdev(msg)
            return
        logger.info(msg)
        msg = datetime.datetime.utcnow().strftime('%m-%d-%Y %H:%M:%S ') + msg
        (output, error, pid) = self.runSubprocess(['slackalert', msg])

    # send debug message to slack
    def slackdebug(self, msg):
        if self.simulate:
            self.slackdev(msg)
            return
        logger.debug(msg)
        msg = datetime.datetime.utcnow().strftime('%m-%d-%Y %H:%M:%S ') + msg
        (output, error, pid) = self.runSubprocess(['slackdebug', msg])

    # send dev message to slack
    def slackdev(self, msg):
        logger.debug(msg)
        msg = datetime.datetime.utcnow().strftime('%m-%d-%Y %H:%M:%S ') + msg
        (output, error, pid) = self.runSubprocess(['slackdev', msg])

    # send preview of fits image to Slack
    def slackpreview(self, fits):
        if self.simulate:
            self.slackdev('Placeholder for image (%s).' % fits)
            return
        (output, error, pid) = self.runSubprocess(
            ['stiffy', fits, 'image.tif'])
        (output, error, pid) = self.runSubprocess(
            ['convert', '-resize', '50%', '-normalize', '-quality', '75', 'image.tif', 'image.jpg'])
        (output, error, pid) = self.runSubprocess(
            ['slackpreview', 'image.jpg', fits])
        (output, error, pid) = self.runSubprocess(
            ['rm', 'image.jpg', 'image.tif'])

    # send preview of fits image to Slack
    def slackimage(self, image):
        if self.simulate:
            self.slackdev('Placeholder for image (%s).' % image)
            return
        (output, error, pid) = self.runSubprocess(
            ['convert', '-resize', '100%', '-normalize', '-quality', '75', image, 'image.jpg'])
        (output, error, pid) = self.runSubprocess(
            ['slackpreview', 'image.jpg', image])
        (output, error, pid) = self.runSubprocess(
            ['rm', 'image.jpg'])

    def squeezeit(self):
        logger.info('Closing the observatory...')
        (output, error, pid) = self.runSubprocess(['squeezeit'], self.simulate)
        # sys.exit(1)

    def crackit(self):
        logger.info('Opening the observatory...')
        (output, error, pid) = self.runSubprocess(
            ['tx', 'lock', 'clear'], self.simulate)
        (output, error, pid) = self.runSubprocess(
            ['tx', 'lock', 'user=mcnowinski', 'email=mcnowinski@gmail.com', 'phone=7032869140', 'comment=chultun'], self.simulate)
        (output, error, pid) = self.runSubprocess(
            ['tin', 'interrupt'], self.simulate)
        (output, error, pid) = self.runSubprocess(
            ['openup', 'nocloud'], self.simulate)
        (output, error, pid) = self.runSubprocess(
            ['keepopen', 'kill', 'debug', 'maxtime=36000', 'slit'], self.simulate, False)
        (output, error, pid) = self.runSubprocess(
            ['tx', 'track', 'on'], self.simulate)
        (output, error, pid) = self.runSubprocess(['tx', 'slit'])
        match = re.search('slit=open', output)
        if match:
            logger.info('Observatory is open.')
            self.is_cracked = True
        else:
            logger.error('Observatory could not be opened.')
        # sys.exit(1)

    def setFocus(self, position):
        if not self.simulate:
            (output, error, pid) = self.runSubprocess(
                ['tx', 'focus', 'pos=%d' % position])

            # done focus pos=4854
            if not re.search('done focus pos\\=([0-9]+)', output):
                logger.error('Could not set focus to %d.' % position)
                return False

        logger.info('Focus position is %s.' % position)
        return True

    def getFocus(self):
        focus_position = ''
        (output, error, pid) = self.runSubprocess(['tx', 'focus'])
        # done focus pos=4854
        match = re.search('pos=(\S+)', output)
        if(match):
            focus_position = match.group(1)
            logger.info('Focus position is %s.' % focus_position)
        else:
            logger.error('Focus command failed (%s).' % output)
        return focus_position

    # check slit
    # if the slit is closed, alert the observer via Slack
    def checkSlit(self):
        (output, error, pid) = self.runSubprocess(['tx', 'slit'])
        match = re.search('slit=open', output)
        if match:
            logger.debug('Slit is open.')
        elif self.simulate:
            logger.debug('Slit is (not really) open.')
        else:
            logger.error(
                'Slit has closed unexpectedly. Trying to re-open it...')
            self.slackalert(
                'Slit has closed unexpectedly. Trying to re-open it...')
            # try to re-open it
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

    # check sun altitude
    # if the sun is too high, squeeze it
    def checkSun(self, wait=False):
        (output, error, pid) = self.runSubprocess(['sun'])
        match = re.search('alt=([\\-\\+\\.0-9]+)', output)
        if match:
            alt = float(match.group(1))
            logger.debug('Sun altitude is %s deg.' % alt)
            if alt > self.max_sun_alt:
                logger.warn('Sun is too high (%s > %s deg).' %
                            (alt, self.max_sun_alt))
                # if wait is True, wait until sun rises
                if wait == True:
                    while alt > self.max_sun_alt:
                        (output, error, pid) = self.runSubprocess(['sun'])
                        match = re.search('alt=([\\-\\+\\.0-9]+)', output)
                        if match:
                            alt = float(match.group(1))
                            self.slackdebug('Sun is too high (%s > %s deg).' %
                                            (alt, self.max_sun_alt))
                        else:
                            logger.error(
                                'Error. Could not determine the current altitude of the sun (%s).' % output)
                            self.slackalert(
                                'Error. Could not determine the current altitude of the sun.')
                            self.done()
                        time.sleep(30)
                else:  # if not wait, we are done!
                    if not self.simulate:
                        self.done()
        else:
            logger.error(
                'Error. Could not determine the current altitude of the sun (%s).' % output)
            self.slackalert(
                'Error. Could not determine the current altitude of the sun.')
            self.done()

    # check altitude of the telescope
    def checkAlt(self):
        (output, error, pid) = self.runSubprocess(['tx', 'where'])
        match = re.search('alt=([\\-\\+\\.0-9]+)', output)
        if match:
            alt = float(match.group(1))
            logger.debug('Telescope altitude is %s deg.' % alt)
            if alt < self.min_alt:
                logger.error(
                    'Telescope altitude is too low (%s < %s deg).' % (alt, self.min_alt))
                self.slackalert(
                    'Telescope altitude is too low (%s < %s deg).' % (alt, self.min_alt))
                if not self.simulate:
                    self.done()
        else:
            logger.error(
                'Error. Could not determine the current altitude of the telescope (%s).' % output)
            self.slackalert(
                'Error. Could not determine the current altitude of the telescope.')
            self.done()

    # check clouds
    # if its too cloudy, wait it out...
    def checkClouds(self):
        (output, error, pid) = self.runSubprocess(['tx', 'taux'])
        match = re.search('cloud=([\\-\\.0-9]+)', output)
        if match:
            clouds = float(match.group(1))
            logger.debug('Cloud cover is %d%%.' % int(clouds*100))
            if clouds >= self.max_clouds_slit:
                logger.error(
                    'Too many clouds (%d%%). Aborting image sequence...' % int(clouds*100))
                if not self.simulate:
                    self.done()
            while clouds >= self.max_clouds_image:
                self.slackalert(
                    'Too many clouds (%d%%). Pausing image sequence...' % int(clouds*100))
                self.checkSun()
                self.checkSlit()
                self.checkAlt()
                if self.simulate:
                    break
                time.sleep(30)
                match = re.search('cloud=([\\-\\.0-9]+)', output)
                if match:
                    clouds = float(match.group(1))
                    if clouds >= self.max_clouds_slit:
                        logger.error(
                            'Too many clouds (%d%%). Aborting image sequence...' % int(clouds*100))
                        if not self.simulate:
                            self.done()
                    logger.debug('Cloud cover is %d%%.' % int(clouds*100))
                else:
                    logger.error('Cloud command failed (%s).' % output)
        return True

    def pinpoint(self, observation, point=True):
        name = observation.target.name
        ra = observation.target.getRa().replace(' ', ':')
        dec = observation.target.getDec().replace(' ', ':')

        self.slackdebug(
            'Pointing telescope to %s (RA=%s, DEC=%s)...' % (name, ra, dec))

        # turn tracking on (just in case)
        (output, error, pid) = self.runSubprocess(
            ['tx', 'track', 'on'], self.simulate)

        # point the telescope
        start_point = datetime.datetime.utcnow()
        if point == True:  # point and refine
            (output, error, pid) = self.runSubprocess(
                ['pinpoint', '%s' % ra, '%s' % dec], self.simulate)
        else:  # just refine
            (output, error, pid) = self.runSubprocess(
                [python_path, pinpoint_path, '%s' % ra, '%s' % dec], self.simulate)
        end_point = datetime.datetime.utcnow()

        # calculate pointing time in seconds
        dt_point = (end_point-start_point).total_seconds()
        logger.debug(
            'Pinpointing telescope required %d seconds to complete.' % dt_point)

        # if refining takes more than max_pinpoint_time_s, send an alert to Slack
        if point == False and dt_point > max_pinpoint_time_s:
            self.slackalert(
                'Warning! Pinpointing telescope required %d seconds to complete. Check clouds, tracking, etc.' % dt_point)

        # check the current telescope position
        (output, error, pid) = self.runSubprocess(['tx', 'where'])

    def pinpointier(self, observation, point=True):
        # MODIFY THESE FIELDS AS NEEDED!
        base_path = '/tmp/'+datetime.datetime.now().strftime("%Y%m%d.%H%M%S%f.pinpoint.")
        # path to astrometry.net solve_field executable
        # astrometry parameters
        downsample = 2
        # bin=1, 0.75 arsec/pixel
        scale_low = 0.55
        scale_high = 2.00
        radius = 30.0  # up this to 30 deg, just in case scope is *way* off
        cpu_limit = 30
        # offset limits (deg)
        max_ra_offset = 30.0
        max_dec_offset = 30.0
        min_ra_offset = 0.05
        min_dec_offset = 0.05
        # how many pointing iterations to allow?
        max_tries = 5
        # image command parameters
        time = 10
        bin = 2
        fits_fname = base_path+'pointing.fits'

        ra_target = Angle(observation.target.getRa().replace(
            ' ', ':'), unit=u.hour).degree
        dec_target = Angle(observation.target.getDec().replace(
            ' ', ':'), unit=u.deg).degree

        # turn tracking on (just in case)
        (output, error, pid) = self.runSubprocess(
            ['tx', 'track', 'on'], self.simulate)

        # make sure dome is properly positioned
        (output, error, pid) = self.runSubprocess(
            ['tx', 'dome', 'center'], self.simulate)

        if point:
            logger.info('Pointing to RA=%s, DEC=%s.' % (observation.target.getRa().replace(
                ' ', ':'),  observation.target.getDec().replace(
                ' ', ':')))
            self.slackdebug('Pointing to RA=%s, DEC=%s.' % (observation.target.getRa().replace(
                ' ', ':'),  observation.target.getDec().replace(
                ' ', ':')))
            (output, error, pid) = self.runSubprocess(
                ['tx', 'point', 'ra=%s' % observation.target.getRa().replace(
                    ' ', ':'), 'dec=%s' % observation.target.getDec().replace(
                    ' ', ':')])

        current_filter = 'clear'
        # ensure filter is clear!
        # get current filter setting
        (output, error, pid) = self.runSubprocess(['pfilter'])
        match = re.search('([a-zA-Z0-1\\-]+)', output)
        if match:
            current_filter = match.group(1)
            logger.debug('Current filter is %s.' % current_filter)
            # set to clear (temporarily)
            logger.debug('Changing filter setting to clear (temporarily).')
            (output, error, pid) = self.runSubprocess(['pfilter', 'clear'])
        else:
            logger.error('Unrecognized filter (%s).' % output)

        ra_offset = 5.0
        dec_offset = 5.0
        iteration = 0
        while((abs(ra_offset) > min_ra_offset or abs(dec_offset) > min_dec_offset) and iteration < max_tries):
            iteration += 1

            logger.debug('Performing adjustment #%d...' % (iteration))

            # get pointing image
            (output, error, pid) = self.runSubprocess(
                ['image', 'time=%f' % time, 'bin=%d' % bin, 'outfile=%s' % fits_fname])

            if not os.path.isfile(fits_fname):
                logger.error('File (%s) not found.' % fits_fname)
                return

            self.slackdebug('Got pinpoint image.')
            self.slackpreview(fits_fname)

            # get FITS header, pull RA and DEC for cueing the plate solving
            if(ra_target == None or dec_target == None):
                header = getheader(fits_fname)
                try:
                    ra_target = header['RA']
                    # create an Angle object
                    ra_target = coord.Angle(ra_target, unit=u.hour).degree
                    dec_target = header['DEC']
                    dec_target = coord.Angle(dec_target, unit=u.deg).degree
                except KeyError:
                    logger.error(
                        "RA/DEC not found in input FITS header (%s)." % fits_fname)
                    return

            # plate solve this image, using RA/DEC from FITS header
            (output, error, pid) = self.runSubprocess([solve_field_path, '--no-verify', '--overwrite', '--no-remove-lines', '--downsample', '%d' % downsample, '--scale-units', 'arcsecperpix', '--no-plots',
                                                       '--scale-low',  '%f' % scale_low, '--scale-high',  '%f' % scale_high, '--ra',  '%s' % ra_target, '--dec', '%s' % dec_target, '--radius',  '%f' % radius, '--cpulimit', '%d' % cpu_limit, fits_fname])
            dumper.debug(output)

            # remove astrometry.net temporary files
            try:
                os.remove(fits_fname)
                os.remove(base_path+'pointing-indx.xyls')
                os.remove(base_path+'pointing.axy')
                os.remove(base_path+'pointing.corr')
                os.remove(base_path+'pointing.match')
                os.remove(base_path+'pointing.rdls')
                os.remove(base_path+'pointing.solved')
                os.remove(base_path+'pointing.wcs')
                os.remove(base_path+'pointing.new')
            except:
                continue

            # look for field center in solve-field output
            match = re.search(
                'Field center\: \(RA,Dec\) \= \(([0-9\-\.\s]+)\,([0-9\-\.\s]+)\) deg\.', output)
            if match:
                RA_image = match.group(1).strip()
                DEC_image = match.group(2).strip()
            else:
                logger.error(
                    "Field center RA/DEC not found in solve-field output!")
                return

            ra_offset = float(ra_target)-float(RA_image)
            if ra_offset > 350:
                ra_offset -= 360.0
            dec_offset = float(dec_target)-float(DEC_image)

            if(abs(ra_offset) <= max_ra_offset and abs(dec_offset) <= max_dec_offset):
                (output, error, pid) = self.runSubprocess(
                    ['tx', 'offset', 'ra=%f' % ra_offset, 'dec=%f' % dec_offset])
                logger.debug("...complete (dRA=%f deg, dDEC=%f deg)." %
                             (ra_offset, dec_offset))
                self.slackdebug("Telescope offset complete (dRA=%f deg, dDEC=%f deg)." % (
                    ra_offset, dec_offset))
            else:
                logger.error("Calculated offsets too large (tx offset ra=%f dec=%f)! Pinpoint aborted." % (
                    ra_offset, dec_offset))
                return

        if(iteration < max_tries):
            logger.info('BAM! Your target has been pinpoint-ed!')
            self.slackdebug('Your target has been pinpoint-ed!')
            return

        logger.error(
            'Exceeded maximum number of adjustments (%d).' % max_tries)
        self.slackalert(
            'Exceeded maximum number of adjustments (%d).' % max_tries)
        return

    def getImage(self, observation, end_time=None):
        if end_time == None:  # no end time, just repeat stack(s) by count
            # make double sure that this is not a continuous observation!
            if observation.sequence.repeat == Sequence.CONTINUOUS:
                logger.error(
                    'Continuous observing sequence received without an end time. Aborting.')
                return
            # how many times does this imaging sequence get repeated?
            for repeat in range(0, observation.sequence.repeat):
                self.getImageStacks(observation)
        else:  # end time specified, this must be a continuous observation
            # make double sure that this is a continuous observation!
            if observation.sequence.repeat != Sequence.CONTINUOUS:
                logger.error(
                    'Non-continuous observing sequence received with an end time. Aborting.')
                return
            logger.debug(
                'Repeating background observation until %s...' % end_time.iso)
            while Time.now() < end_time:  # how many times does this imaging sequence get repeated?
                self.getImageStacks(observation)

    def image(self, exposure, filter, binning, path, filename):
        # make sure path exists!
        pathlib2.Path('%s' % (path)).mkdir(parents=True, exist_ok=True)
        # get the image!
        if filter == 'dark':
            # use h-alpha filter to reduce any ambient light
            (output, error, pid) = self.runSubprocess(
                ['pfilter', 'h-alpha'])
            (output, error, pid) = self.runSubprocess(
                ['image', 'dark', 'time=%f' % exposure, 'bin=%d' % binning, 'outfile=%s/%s' % (path, filename)])
        else:
            (output, error, pid) = self.runSubprocess(
                ['pfilter', '%s' % filter])
            (output, error, pid) = self.runSubprocess(
                ['image', 'time=%f' % exposure, 'bin=%d' % binning, 'outfile=%s/%s' % (path, filename)])

    def getImageStacks(self, observation, end_time=None, doChecks=True):
        for stack in observation.sequence.stacks:  # iterate thru the image sets
            # how many times does this stack get repeated?
            for count in range(0, stack.count):
                if doChecks:
                    # check sun, clouds, slit, etc.
                    self.checkSun()
                    if self.is_cracked:
                        self.checkSlit()
                    self.checkAlt()
                    self.checkClouds()

                name = observation.target.name.strip()
                name = name.replace(' ', '_')
                name = name.replace('(', '')
                name = name.replace(')', '')
                name = name.replace("'", '')
                filter = stack.filter
                exposure = stack.exposure
                binning = stack.binning
                do_pinpoint = stack.do_pinpoint
                user = observation.observer
                # get image
                # M42_clear_100s_bin1_180210_043012_seo_epjmm_003_RAW.fits
                # fits = '%s_%s_%dsec_bin%d_%s_%s_num%d_seo.fits' % (
                #    name, filter, exposure, binning, user, datetime.datetime.utcnow().strftime('%Y%b%d_%Hh%Mm%Ss'), count)
                fits = '%s_%s_%.2fs_bin%d_%s_seo_%s_%03d_RAW.fits' % (
                    name, filter, exposure, binning, datetime.datetime.utcnow().strftime('%y%m%d_%H%M%S'), user, count)
                fits = fits.replace(' ', '_')
                fits = fits.replace('(', '')
                fits = fits.replace(')', '')
                fits = fits.replace("'", '')
                self.slackdebug('Taking image (%s)...' % (fits))
                if self.simulate:
                    time.sleep(exposure)
                    error = 0
                else:
                    self.image(exposure, filter, binning,
                               image_path+'/'+user+'/'+name, fits)
#                    # make sure path exists!
#                    pathlib2.Path('%s/%s/%s' % (image_path, user, name)
#                                  ).mkdir(parents=True, exist_ok=True)
#                    # get the image!
#                    if filter == 'dark':
#                        # use h-alpha filter to reduce any ambient light
#                        (output, error, pid) = self.runSubprocess(
#                            ['pfilter', 'h-alpha'])
#                        (output, error, pid) = self.runSubprocess(
#                            ['image', 'dark', 'time=%f' % exposure, 'bin=%d' % binning, 'outfile=%s/%s/%s/%s' % (image_path, user, name, fits)])
#                    else:
#                        (output, error, pid) = self.runSubprocess(
#                            ['pfilter', '%s' % filter])
#                        (output, error, pid) = self.runSubprocess(
#                            ['image', 'time=%f' % exposure, 'bin=%d' % binning, 'outfile=%s/%s/%s/%s' % (image_path, user, name, fits)])
                # if not error:
                self.slackdebug('Got image (%s/%s/%s/%s).' %
                                (image_path, user, name, fits))
                if not self.simulate:
                    self.slackpreview('%s/%s/%s/%s' %
                                      (image_path, user, name, fits))
                    # IMAGE_PATHNAME=$STARS_IMAGE_PATH/`date -u +"%Y"`/`date -u +"%Y-%m-%d"`/${NAME}
                    #(ssh -q -i $STARS_PRIVATE_KEY_PATH $STARS_USERNAME@$STARS_SERVER "mkdir -p $IMAGE_PATHNAME"; scp -q -i $STARS_PRIVATE_KEY_PATH $IMAGE_FILENAME $STARS_USERNAME@$STARS_SERVER:$IMAGE_PATHNAME/$IMAGE_FILENAME) &
                    #(output, error, pid) = runSubprocess(['tostars','%s'%name.replace(' ', '_').replace('(', '').replace(')', ''),'%s'%fits])
                else:
                    self.slackdebug(
                        'Error. Image command failed (%s/%s/%s/%s).' % (image_path, user, name, fits))
                if do_pinpoint:
                    self.pinpointier(observation, False)

    # ssh -q -i /home/mcnowinski/.ssh/id_rsa_stars dmcginnis427@stars.uchicago.edu "mkdir -p /data/images/StoneEdge/0.5meter/2018/2018-05-26/ exit"
    # scp -q -r -i /home/mcnowinski/.ssh/id_rsa_stars * dmcginnis427@stars.uchicago.edu:/data/images/StoneEdge/0.5meter/2018/2018-05-26/
    def toStars(self):
        # are there any .fits images to send?
        fits = glob2.glob(image_path+'/**/*.fits')
        if(len(fits) <= 0):
            logger.info('No fits files found in %s.' % image_path)
            return
        # we are going to put these in a folder corresponding to the datetime this command was run!
        # a bit different from how this usually works...
        dest_path = '%s%s/%s' % (stars_image_path, datetime.datetime.utcnow().strftime(
            '%Y'), datetime.datetime.utcnow().strftime('%Y-%m-%d'))
        # off they go!
        # create new directory if required
        #(output, error, pid) = runSubprocess(
        #    ['ssh', '-q', '-i', %s*.fits' % image_path, '%s' % dest_path])
        # if error == '':
        #    logger.info('Successfully uploaded %d image(s) to stars (%s).' % (
        #        len(fits), dest_path))
        #    # move images to archive
        #    files = glob.iglob(image_path+'*.fits')
        #    for file in files:
        #        if os.path.isfile(file):
        #            shutil.move(file, image_archive_path)
        # else:
        #    send_message('Error. Image uploaded failed!')
#
# an astronomical observation
#


class Observation():

    sequence = None  # image sequence
    target = None  # target
    observatory = None  # observatory
    observer = None  # who's observing?
    min_obs_alt = None  # min alt to start observations in deg

    # used by the Scheduler and others
    obs_start_time = None  # when will target best be observable?
    min_obs_time = None  # when is the target first observable?
    max_obs_time = None  # when is the target last observable?
    max_alt_time = None  # when is the
    active = True  # is this observation (still) active?
    id = -1  # id

    # init
    def __init__(self, observatory, target, sequence, min_obs_alt, observer):
        self.observatory = observatory
        self.target = target
        self.sequence = sequence
        self.min_obs_alt = min_obs_alt
        self.observer = observer
        # self.getTimes()

    def toString(self):
        return '%s\n%s\n%smin_alt=%f deg\nobs_time=%s\nid=%d\nactive=%d\nmin_obs_time=%s\nmax_obs_time=%s\nmax_alt_time=%s\nuser=%s' % (self.observatory.toString(), self.target.toString(), self.sequence.toString(), self.min_obs_alt, self.obs_start_time, self.id, self.active, self.min_obs_time, self.max_obs_time, self.max_alt_time, self.observer)

    # for this observation, get min/max observable times and max alt time
    def getTimes(self):
        # temp var to hold obs info
        obs = {'time': [], 'id': []}

        # init observatory location
        observatory_location = EarthLocation(
            lat=self.observatory.latitude*u.deg, lon=self.observatory.longitude*u.deg, height=self.observatory.altitude*u.m)
        # get next sunrise and nearest sunset times
        observatory_location_obsplan = Observer(longitude=self.observatory.longitude*u.deg, latitude=self.observatory.latitude *
                                                u.deg, elevation=self.observatory.altitude*u.m, name=self.observatory.code, timezone=self.observatory.timezone)
        sunset_time = observatory_location_obsplan.twilight_evening_nautical(
            Time.now(), which="nearest")
        sunrise_time = observatory_location_obsplan.twilight_morning_nautical(
            Time.now(), which="next")
        logger.debug('The nearest sunset is %s. The next sunrise is %s.' %
                     (sunset_time.iso, sunrise_time.iso))

        # build alt-az coordinate frame for observatory over next ? hours (e.g., nearest sunset to next sunrise)
        # start time is sunset or current time, if later...
        now = Time.now()
        if (now > sunset_time):
            obs_time = Time.now()
        else:
            obs_time = sunset_time
        delta_obs_time = np.linspace(
            0, (sunrise_time-obs_time).sec/3600., 1000)*u.hour
        # array of times between sunset and sunrise
        times = obs_time + delta_obs_time
        # celestial frame for this observatory over times
        frame = AltAz(obstime=times, location=observatory_location)

        # build target altaz relative to observatory
        target_ra = self.target.getRa()
        target_dec = self.target.getDec()
        input_coordinates = target_ra + " " + target_dec
        try:
            target_coordinates = SkyCoord(
                input_coordinates, unit=(u.hourangle, u.deg))
        except:
            pass
        target_altaz = target_coordinates.transform_to(frame)

        # when is target highest *and* above minimum altitude?
        # when is it above min_obs_alt?
        valid_alt_times = times[np.where(
            target_altaz.alt >= self.min_obs_alt*u.degree)]
        # when does the max alt occur?
        if len(valid_alt_times) > 0:
            self.min_obs_time = Time(
                np.min(times[np.where(target_altaz.alt > self.min_obs_alt*u.degree)]))
            self.max_obs_time = Time(
                np.max(times[np.where(target_altaz.alt > self.min_obs_alt*u.degree)]))
            self.max_alt_time = Time(
                times[np.argmax(target_altaz.alt)])
        else:
            logger.error('Target (%s) is not observable.' %
                         self.target.getName())

#
# the astronomical target
#


class Target():

    # init
    def __init__(self, name, ra, dec):
        self.name = name
        self.ra = ra  # hour:min:sec
        self.dec = dec  # deg:min:sec

    # init with name and type only
    @classmethod
    def from_name(cls, keyword, observatory, type):
        objects = Target.findObjects(keyword, observatory, type)
        if len(objects) == 0:
            logger.error('Could not find matching object for %s.' % keyword)
            sys.exit(1)
        else:
            if len(objects) > 1:
                logger.warn('Found multiple matching objects for %s. Using first object (%s).' % (
                    name, objects[0]['name']))
        target = cls(objects[0]['name'], objects[0]['ra'], objects[0]['dec'])
        return target

    # name
    def getName(self):
        return self.name

    def setName(self, name):
        self.name = name

    # ra = right ascension
    def getRa(self):
        return self.ra

    def setRa(self, ra):
        self.ra = ra

    #dec = declination
    def getDec(self):
        return self.dec

    def setDec(self, dec):
        self.dec = dec

    def toString(self):
        return 'target: name=%s, ra=%s, dec=%s' % (self.name, self.ra, self.dec)

    @staticmethod
    def findObjects(keyword, observatory, type):
        type = type.lower()
        if (type == 'asteroid' or type == 'planet' or type == 'solar system'):
            return Target.findSolarSystemObjects(keyword, observatory)
        elif (type == 'star' or type == 'celestial' or type == 'galaxy'):
            return Target.findCelestialObjects(keyword)
        else:
            logger.error("Unknown type (%s) in Target.findObjects." % type)
            return []

    @staticmethod
    def findCelestialObjects(keyword):
        results = Simbad.query_object(keyword)
        if results == None:
            return []
        objects = []
        for result in results:
            objects.append({'type': 'Celestial', 'id': result['MAIN_ID'], 'name': result['MAIN_ID'].replace(' ', ''),
                            'ra': result['RA'], 'dec': result['DEC']})
        return objects

    # search solar system small bodies using JPL HORIZONS
    @staticmethod
    def findSolarSystemObjects(keyword, observatory):
        # ch constants
        # max airmass
        max_airmass = 2.0  # 30 deg elevation
        objects = []
        # list of matches
        object_names = []
        # set to * to make the searches wider by default
        suffix = ''
        # two passes, one for major (and maybe small) and one for (only) small bodies
        lookups = [keyword + suffix, keyword + suffix + ';']
        for repeat in range(0, 2):
            # user JPL Horizons batch to find matches
            f = urllib2.urlopen('https://ssd.jpl.nasa.gov/horizons_batch.cgi?batch=l&COMMAND="%s"' %
                                urllib2.quote(lookups[repeat].upper()))
            output = f.read()  # the whole enchilada
            #print output
            lines = output.splitlines()  # line by line
            # no matches? go home
            if re.search('No matches found', output):
                logger.debug('No matches found in JPL Horizons for %s.' %
                             lookups[repeat].upper())
            elif re.search('Target body name:', output):
                logger.debug('Single match found in JPL Horizons for %s.' %
                             lookups[repeat].upper().replace(suffix, ''))
                # just one match?
                # if major body search (repeat = 0), ignore small body results
                # if major body search, grab integer id
                if repeat == 0:
                    if re.search('Small-body perts:', output):
                        continue
                    match = re.search(
                        'Target body name:\\s[a-zA-Z]+\\s\\((\\d+)\\)', output)
                    if match:
                        object_names.append(match.group(1))
                    else:
                        logger.error('Error. Could not parse id for single match major body (%s).' %
                                     lookups[repeat].upper().replace(suffix, ''))
                else:
                    # user search term is unique, so use it!
                    object_names.append(
                        lookups[repeat].upper().replace(suffix, ''))
            elif repeat == 1 and re.search('Matching small-bodies', output):
                logger.info('Multiple small bodies found in JPL Horizons for %s.' %
                            lookups[repeat].upper())
                # Matching small-bodies:
                #
                #    Record #  Epoch-yr  Primary Desig  >MATCH NAME<
                #    --------  --------  -------------  -------------------------
                #          4             (undefined)     Vesta
                #      34366             2000 RP36       Rosavestal
                match_count = 0
                for line in lines:
                    search_string = line.strip()
                    # look for small body list
                    match = re.search('^-?\\d+', search_string)
                    # parse out the small body parameters
                    if match:
                        match_count += 1
                        record_number = line[0:12].strip()
                        epoch_yr = line[12:22].strip()
                        primary_desig = line[22:37].strip()
                        match_name = line[37:len(line)].strip()
                        #print record_number, epoch_yr, primary_desig, match_name
                        # add semicolon for small body lookups
                        object_names.append(record_number + ';')
                # check our parse job
                match = re.search('(\\d+) matches\\.', output)
                if match:
                    if int(match.group(1)) != match_count:
                        logger.error('Multiple JPL small body parsing error!')
                    else:
                        logger.info(
                            'Multiple JPL small body parsing successful!')
            elif repeat == 0 and re.search('Multiple major-bodies', output):
                logger.info('Multiple major bodies found in JPL Horizons for %s.' %
                            lookups[repeat].upper())
                # Multiple major-bodies match string "50*"
                #
                #  ID#      Name                               Designation  IAU/aliases/other
                #  -------  ---------------------------------- -----------  -------------------
                #      501  Io                                              JI
                #      502  Europa                                          JII
                match_count = 0
                for line in lines:
                    search_string = line.strip()
                    # look for major body list
                    match = re.search('^-?\\d+', search_string)
                    # parse out the major body parameters
                    if match:
                        match_count += 1
                        record_number = line[0:9].strip()
                        # negative major bodies are spacecraft,etc. Skip those!
                        if int(record_number) >= 0:
                            name = line[9:45].strip()
                            designation = line[45:57].strip()
                            other = line[57:len(line)].strip()
                            #print record_number, name, designation, other
                            # NO semicolon for major body lookups
                            object_names.append(record_number)
                # check our parse job
                match = re.search('Number of matches =([\\s\\d]+).', output)
                if match:
                    if int(match.group(1)) != match_count:
                        logger.error('Multiple JPL major body parsing error!')
                    else:
                        logger.info(
                            'Multiple JPL major body parsing successful!')
        # get *nearest* sunset and *next* sunrise times
        # still not a big fan of this!
        observatory_location_obsplan = Observer(longitude=observatory.longitude*u.deg, latitude=observatory.latitude *
                                                u.deg, elevation=observatory.altitude*u.m, name=observatory.code, timezone=observatory.timezone)
        start = observatory_location_obsplan.twilight_evening_nautical(
            Time.now(), which="nearest")
        end = observatory_location_obsplan.twilight_morning_nautical(
            Time.now(), which="next")
        logger.debug('The nearest sunset is %s. The next sunrise is %s.' %
                     (start.iso, end.iso))
        logger.info('Found %d solar system match(es) for "%s".' %
                    (len(object_names), keyword))
        count = 0
        for object_name in object_names:
            count += 1
            # get ephemerides for target in JPL Horizons from start to end times
            result = ch.query(object_name.upper(), smallbody=True)
            result.set_epochrange(start.iso, end.iso, '15m')
            result.get_ephemerides(observatory.code)
            # return transit RA/DEC if available times exist
            if result and len(result['EL']):
                imax = np.argmax(result['EL'])
                ra = Angle(float(result['RA'][imax]) *
                           u.deg).to_string(unit=u.hour, sep=':')
                dec = Angle(float(result['DEC'][imax]) *
                            u.deg).to_string(unit=u.degree, sep=':')
                objects.append({'type': 'Solar System', 'id': object_name.upper(
                ), 'name': result['targetname'][0], 'ra': ra, 'dec': dec})
            else:
                logger.debug('The object ('+object_name+') is not observable.')
        return objects

#
# settings for a single set of astronomical images
#


class Stack():

    exposure = 10  # exposure time in seconds
    filter = 'clear'  # filter, e.g., clear, h-alpha, u-band, g-band, r-band, i-band, z-band
    binning = 1  # binning, e.g. 1 or 2
    count = 1  # number of images in this stack
    do_pinpoint = True  # refine pointing in between imaging

    # init
    def __init__(self, exposure, filter, binning, count, do_pinpoint=True):
        self.exposure = exposure
        self.filter = filter
        self.binning = binning
        self.count = count
        self.do_pinpoint = do_pinpoint

    def toString(self):
        return 'image stack: exposure=%f, filter=%s, binning=%d, count=%d' % (self.exposure, self.filter, self.binning, self.count)


#
# sequence of astronomical image stacks
#
class Sequence():

    stacks = []  # list of image stacks
    repeat = None  # number of times to repeat this sequence

    # repeat as much as possible
    CONTINUOUS = -1

    # init
    def __init__(self, stacks, repeat):
        self.stacks = stacks
        self.repeat = repeat

    def addStack(self, stack):
        self.stacks.append(stack)

    def toString(self):
        sequence_string = 'sequence: repeat=%d\n' % (self.repeat)
        for stack in self.stacks:
            sequence_string += '  %s\n' % stack.toString()
        return sequence_string

#
# the observatory
#


class Observatory():

    code = None
    latitude = 0.0  # in decimal degrees
    longitude = 0.0  # in decimal degrees
    altitude = 0.0  # in meters
    timzeone = None

    # init
    def __init__(self, code, latitude, longitude, altitude, timezone):
        self.code = code
        self.latitude = latitude
        self.longitude = longitude
        self.altitude = altitude
        self.timezone = timezone

    def toString(self):
        observatory_string = 'observatory: code=%s, lat=%f, lon=%f, alt=%f' % (
            self.code, self.latitude, self.longitude, self.altitude)
        return observatory_string


class Scheduler():
    # based on Amanda Pagul's schedule function for the SEO queue
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
    observatory = None  # an Observatory
    observations = None  # list of Observations
    last_id = 0  # unique id for each observation

    sunset_time = None  # nearest sunset
    sunrise_time = None  # next sunrise

    # max_sun_alt = -12 #what defines dark? (deg)
    # min_target_alt = 30 #how low can you go? (deg)

    def __init__(self, observatory, observations=None):
        self.observatory = observatory

        # ignore warnings
        warnings.filterwarnings('ignore', category=UserWarning, append=True)
        warnings.filterwarnings('ignore', category=FutureWarning, append=True)

        # uncomment to download the latest for astroplan
        logger.debug('Updating Astroplan IERS Bulletin A...')
        try:
            download_IERS_A()
        except:
            logger.error('Error. Could not update Astroplan IERS Bulletin A.')

        # get *nearest* sunset and *next* sunrise times
        # still not a big fan of this!
        observatory_location_obsplan = Observer(longitude=self.observatory.longitude*u.deg, latitude=self.observatory.latitude *
                                                u.deg, elevation=self.observatory.altitude*u.m, name=self.observatory.code, timezone=self.observatory.timezone)
        self.sunset_time = observatory_location_obsplan.twilight_evening_nautical(
            Time.now(), which="nearest")
        self.sunrise_time = observatory_location_obsplan.twilight_morning_nautical(
            Time.now(), which="next")
        logger.debug('The nearest sunset is %s. The next sunrise is %s.' %
                     (self.sunset_time.iso, self.sunrise_time.iso))

        # add observations if provided
        if observations != None:
            self.addObservations(observations)

    def addObservations(self, observations):
        self.observations = observations
        # assign unique ids to observations
        for observation in self.observations:
            observation.id = self.last_id
            self.last_id += 1

    # mark this observation as complete
    def isDone(self, observation):
        for obs in self.observations:
            if obs.id == observation.id:
                obs.active = False
                return
        logger.error('Matching observation (%d) not found.' % observation.id)

    # grab the next (best) observation from the list
    def whatsNext(self):
        # temp var to hold obs info
        obs = {'time': [], 'id': []}

        # init observatory location
        observatory_location = EarthLocation(
            lat=self.observatory.latitude*u.deg, lon=self.observatory.longitude*u.deg, height=self.observatory.altitude*u.m)

        # build alt-az coordinate frame for observatory over next ? hours (e.g., nearest sunset to next sunrise)
        # start time is sunset or current time, if later...
        now = Time.now()
        if (now > self.sunset_time):
            obs_time = Time.now()
        else:
            obs_time = self.sunset_time
        delta_obs_time = np.linspace(
            0, (self.sunrise_time-obs_time).sec/3600., 1000)*u.hour
        # array of times between sunset and sunrise
        times = obs_time + delta_obs_time
        # celestial frame for this observatory over times
        frame = AltAz(obstime=times, location=observatory_location)

        # loop thru observations, suggest the next best target based on time of max alt.
        for observation in self.observations:

            # skip obs that are complete/inactive
            if observation.active == False:
                logger.debug('Observation (%s) is not active. Skipping...' %
                             observation.target.getName())
                continue

            # skip background obs
            if observation.sequence.repeat == Sequence.CONTINUOUS:
                logger.debug('Observation (%s) is not foreground type. Skipping...' %
                             observation.target.getName())
                continue

            # build target altaz relative to observatory
            target_ra = observation.target.getRa()
            target_dec = observation.target.getDec()
            input_coordinates = target_ra + " " + target_dec
            try:
                target_coordinates = SkyCoord(
                    input_coordinates, unit=(u.hourangle, u.deg))
            except:
                continue
            target_altaz = target_coordinates.transform_to(frame)

            # when is target highest *and* above minimum altitude?
            # when is it above min_obs_alt?
            valid_alt_times = times[np.where(
                target_altaz.alt >= observation.min_obs_alt*u.degree)]
            # when does the max alt occur?
            if len(valid_alt_times) > 0:
                obs['id'].append(observation.id)
                # min time is the selection criteria
                obs['time'].append(times[np.argmax(target_altaz.alt)])
                # set min and max obs times
                observation.min_obs_time = Time(
                    np.min(times[np.where(target_altaz.alt > observation.min_obs_alt*u.degree)]))
                observation.max_obs_time = Time(
                    np.max(times[np.where(target_altaz.alt > observation.min_obs_alt*u.degree)]))

        # get earliest max alt. from valid targets
        if len(obs['time']) > 0:
            # the winner!
            id = np.argmin(obs['time'])
            # set obs start time
            self.observations[obs['id'][id]
                              ].obs_start_time = Time(obs['time'][id])
        else:
            logger.info("No active foreground observations exist.")
            return None

        logger.debug('Next foreground observation is %s.' %
                     self.observations[obs['id'][id]].target.name)
        return self.observations[obs['id'][id]]

    # grab the highest observation from the list
    def whatsHighest(self):
        # temp var to hold obs info
        obs = {'time': [], 'delta_time': [], 'id': []}

        # init observatory location
        observatory_location = EarthLocation(
            lat=self.observatory.latitude*u.deg, lon=self.observatory.longitude*u.deg, height=self.observatory.altitude*u.m)

        # build alt-az coordinate frame for observatory over next ? hours (e.g., nearest sunset to next sunrise)
        # start time is sunset or current time, if later...
        now = Time.now()
        if (now > self.sunset_time):
            obs_time = Time.now()
        else:
            obs_time = self.sunset_time
        delta_obs_time = np.linspace(
            0, (self.sunrise_time-obs_time).sec/3600., 1000)*u.hour
        # array of times between sunset and sunrise
        times = obs_time + delta_obs_time
        # celestial frame for this observatory over times
        frame = AltAz(obstime=times, location=observatory_location)

        # loop thru observations, suggest the next best target based on time of max alt.
        for observation in self.observations:

            # skip obs that are complete/inactive
            if observation.active == False:
                logger.debug('Observation (%s) is not active. Skipping...' %
                             observation.target.getName())
                continue

            # skip background obs
            if observation.sequence.repeat == Sequence.CONTINUOUS:
                logger.debug('Observation (%s) is not foreground type. Skipping...' %
                             observation.target.getName())
                continue

            # build target altaz relative to observatory
            target_ra = observation.target.getRa()
            target_dec = observation.target.getDec()
            input_coordinates = target_ra + " " + target_dec
            try:
                target_coordinates = SkyCoord(
                    input_coordinates, unit=(u.hourangle, u.deg))
            except:
                continue
            target_altaz = target_coordinates.transform_to(frame)

            # when is target highest *and* above minimum altitude?
            # when is it above min_obs_alt?
            valid_alt_times = times[np.where(
                target_altaz.alt >= observation.min_obs_alt*u.degree)]
            # when does the max alt occur?
            if len(valid_alt_times) > 0:
                obs['id'].append(observation.id)
                # min time is the selection criteria
                obs['time'].append(times[np.argmax(target_altaz.alt)])
                delta_time = times[np.argmax(target_altaz.alt)]-Time.now()
                delta_time = abs(delta_time.sec)
                obs['delta_time'].append(delta_time)
                # set min and max obs times
                observation.min_obs_time = Time(
                    np.min(times[np.where(target_altaz.alt > observation.min_obs_alt*u.degree)]))
                observation.max_obs_time = Time(
                    np.max(times[np.where(target_altaz.alt > observation.min_obs_alt*u.degree)]))

        # get earliest max alt. from valid targets
        if len(obs['time']) > 0:
            # the winner!
            id = np.argmin(obs['delta_time'])
            # set obs start time
            self.observations[obs['id'][id]
                              ].obs_start_time = Time(obs['time'][id])
        else:
            logger.info("No active foreground observations exist.")
            return None

        logger.debug('Next foreground observation is %s.' %
                     self.observations[obs['id'][id]].target.name)
        return self.observations[obs['id'][id]]

    # grab the next (best) background observation
    def whatsInBetween(self):
        # temp var to hold obs info
        obs = {'time': [], 'id': []}

        # init observatory location
        observatory_location = EarthLocation(
            lat=self.observatory.latitude*u.deg, lon=self.observatory.longitude*u.deg, height=self.observatory.altitude*u.m)

        # build alt-az coordinate frame for observatory over next ? hours (e.g., nearest sunset to next sunrise)
        # start time is sunset or current time, if later...
        now = Time.now()
        if (now > self.sunset_time):
            obs_time = Time.now()
        else:
            obs_time = self.sunset_time
        delta_obs_time = np.linspace(
            0, (self.sunrise_time-obs_time).sec/3600., 1000)*u.hour
        # array of times between sunset and sunrise
        times = obs_time + delta_obs_time
        # celestial frame for this observatory over times
        frame = AltAz(obstime=times, location=observatory_location)

        # loop thru observations, suggest the next best target based on time of max alt.
        for observation in self.observations:

            # skip obs that are complete/inactive
            if observation.active == False:
                logger.debug('Observation (%s) is not active. Skipping...' %
                             observation.target.getName())
                continue

            # skip foreground obs
            if observation.sequence.repeat != Sequence.CONTINUOUS:
                logger.debug('Observation (%s) is not background type. Skipping...' %
                             observation.target.getName())
                continue

            # build target altaz relative to observatory
            target_ra = observation.target.getRa()
            target_dec = observation.target.getDec()
            input_coordinates = target_ra + " " + target_dec
            try:
                target_coordinates = SkyCoord(
                    input_coordinates, unit=(u.hourangle, u.deg))
            except:
                continue
            target_altaz = target_coordinates.transform_to(frame)

            # when is target above minimum altitude?
            # when is it above min_obs_alt?
            valid_alt_times = times[np.where(
                target_altaz.alt >= observation.min_obs_alt*u.degree)]
            # when does the max alt occur?
            if len(valid_alt_times) > 0:
                # set min and max obs times
                observation.min_obs_time = Time(
                    np.min(times[np.where(target_altaz.alt > observation.min_obs_alt*u.degree)]))
                observation.max_obs_time = Time(
                    np.max(times[np.where(target_altaz.alt > observation.min_obs_alt*u.degree)]))
                observation.max_alt_time = Time(
                    times[np.argmax(target_altaz.alt)])
                obs['id'].append(observation.id)
                # remember the end time too
                obs['time'].append(
                    np.max(times[np.where(target_altaz.alt > observation.min_obs_alt*u.degree)]))

        # get earliest max. time from valid targets
        if len(obs) > 0:
            # the winner!
            id = np.argmin(obs['time'])
            # set obs start time
            self.observations[obs['id'][id]].obs_start_time = Time.now()
        else:
            logger.info("No active background observations exist.")
            return None

        logger.debug('Next background observation is %s.' %
                     self.observations[obs['id'][id]].target.name)
        return self.observations[obs['id'][id]]

    def findObjects(self, keyword):
        results = Simbad.query_object(keyword)
        if results == None:
            return []
        objects = []
        for result in results:
            objects.append({'type': 'Celestial', 'id': result['MAIN_ID'], 'name': result['MAIN_ID'].replace(' ', ''),
                            'RA': result['RA'], 'DEC': result['DEC']})
        return objects

    # search solar system small bodies using JPL HORIZONS
    def findSolarSystemObjects(self, keyword):
        # ch constants
        # max airmass
        max_airmass = 2.0  # 30 deg elevation
        objects = []
        # list of matches
        object_names = []
        # set to * to make the searches wider by default
        suffix = ''
        # two passes, one for major (and maybe small) and one for (only) small bodies
        lookups = [keyword + suffix, keyword + suffix + ';']
        for repeat in range(0, 2):
            # user JPL Horizons batch to find matches
            f = urllib2.urlopen('https://ssd.jpl.nasa.gov/horizons_batch.cgi?batch=l&COMMAND="%s"' %
                                urllib2.quote(lookups[repeat].upper()))
            output = f.read()  # the whole enchilada
            #print output
            lines = output.splitlines()  # line by line
            # no matches? go home
            if re.search('No matches found', output):
                logger.info('No matches found in JPL Horizons for %s.' %
                            lookups[repeat].upper())
            elif re.search('Target body name:', output):
                logger.info('Single match found in JPL Horizons for %s.' %
                            lookups[repeat].upper().replace(suffix, ''))
                # just one match?
                # if major body search (repeat = 0), ignore small body results
                # if major body search, grab integer id
                if repeat == 0:
                    if re.search('Small-body perts:', output):
                        continue
                    match = re.search(
                        'Target body name:\\s[a-zA-Z]+\\s\\((\\d+)\\)', output)
                    if match:
                        object_names.append(match.group(1))
                    else:
                        logger.error('Error. Could not parse id for single match major body (%s).' %
                                     lookups[repeat].upper().replace(suffix, ''))
                else:
                    # user search term is unique, so use it!
                    object_names.append(
                        lookups[repeat].upper().replace(suffix, ''))
            elif repeat == 1 and re.search('Matching small-bodies', output):
                logger.info('Multiple small bodies found in JPL Horizons for %s.' %
                            lookups[repeat].upper())
                # Matching small-bodies:
                #
                #    Record #  Epoch-yr  Primary Desig  >MATCH NAME<
                #    --------  --------  -------------  -------------------------
                #          4             (undefined)     Vesta
                #      34366             2000 RP36       Rosavestal
                match_count = 0
                for line in lines:
                    search_string = line.strip()
                    # look for small body list
                    match = re.search('^-?\\d+', search_string)
                    # parse out the small body parameters
                    if match:
                        match_count += 1
                        record_number = line[0:12].strip()
                        epoch_yr = line[12:22].strip()
                        primary_desig = line[22:37].strip()
                        match_name = line[37:len(line)].strip()
                        #print record_number, epoch_yr, primary_desig, match_name
                        # add semicolon for small body lookups
                        object_names.append(record_number + ';')
                # check our parse job
                match = re.search('(\\d+) matches\\.', output)
                if match:
                    if int(match.group(1)) != match_count:
                        logger.error('Multiple JPL small body parsing error!')
                    else:
                        logger.info(
                            'Multiple JPL small body parsing successful!')
            elif repeat == 0 and re.search('Multiple major-bodies', output):
                logger.info('Multiple major bodies found in JPL Horizons for %s.' %
                            lookups[repeat].upper())
                # Multiple major-bodies match string "50*"
                #
                #  ID#      Name                               Designation  IAU/aliases/other
                #  -------  ---------------------------------- -----------  -------------------
                #      501  Io                                              JI
                #      502  Europa                                          JII
                match_count = 0
                for line in lines:
                    search_string = line.strip()
                    # look for major body list
                    match = re.search('^-?\\d+', search_string)
                    # parse out the major body parameters
                    if match:
                        match_count += 1
                        record_number = line[0:9].strip()
                        # negative major bodies are spacecraft,etc. Skip those!
                        if int(record_number) >= 0:
                            name = line[9:45].strip()
                            designation = line[45:57].strip()
                            other = line[57:len(line)].strip()
                            #print record_number, name, designation, other
                            # NO semicolon for major body lookups
                            object_names.append(record_number)
                # check our parse job
                match = re.search('Number of matches =([\\s\\d]+).', output)
                if match:
                    if int(match.group(1)) != match_count:
                        logger.error('Multiple JPL major body parsing error!')
                    else:
                        logger.info(
                            'Multiple JPL major body parsing successful!')
        #start = datetime.datetime.utcnow()
        #end = start+datetime.timedelta(hours=12)
        # get *nearest* sunset and *next* sunrise times
        # still not a big fan of this!
        observatory_location_obsplan = Observer(longitude=self.observatory.longitude*u.deg, latitude=self.observatory.latitude *
                                                u.deg, elevation=self.observatory.altitude*u.m, name=self.observatory.code, timezone=self.observatory.timezone)
        start = observatory_location_obsplan.twilight_evening_nautical(
            Time.now(), which="nearest")
        end = observatory_location_obsplan.twilight_morning_nautical(
            Time.now(), which="next")
        logger.debug('The nearest sunset is %s. The next sunrise is %s.' %
                     (start.iso, end.iso))
        logger.info('Found %d solar system match(es) for "%s".' %
                    (len(object_names), keyword))
        count = 0
        for object_name in object_names:
            count += 1
            # get ephemerides for target in JPL Horizons from start to end times
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
            result = ch.query(object_name.upper(), smallbody=True)
            result.set_epochrange(start.iso, end.iso, '15m')
            # result.get_ephemerides(self.observatory.code,
            #                   skip_daylight=True, airmass_lessthan=max_airmass)
            result.get_ephemerides(self.observatory.code)
            # return transit RA/DEC if available times exist
            if result and len(result['EL']):
                imax = np.argmax(result['EL'])
                ra = Angle(float(result['RA'][imax]) *
                           u.deg).to_string(unit=u.hour, sep=':')
                dec = Angle(float(result['DEC'][imax]) *
                            u.deg).to_string(unit=u.degree, sep=':')
                objects.append({'type': 'Solar System', 'id': object_name.upper(
                ), 'name': result['targetname'][0], 'RA': ra, 'DEC': dec})
            else:
                logger.debug('The object ('+object_name+') is not observable.')
        return objects
