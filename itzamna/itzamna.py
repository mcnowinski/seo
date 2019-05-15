# Web and Real-Time Messaging interface to Slack for use with the SEO telescope, a.k.a. itzamna

from slackclient import SlackClient
import traceback
import time
import random
import string
import math
import os
import sys
import datetime
import json
import re
import urllib2
from StringIO import StringIO
from zipfile import ZipFile
import requests
import subprocess
from astropy.coordinates import SkyCoord, EarthLocation, AltAz, Angle, get_sun
import astropy.units as u
from astropy.time import Time
# use custom callhorizons
import ch  # callhorizons module edited by me
import numpy
import ephem  # calculate satellite (natural and artifical) ephemerides
from astroquery.simbad import Simbad
# for plots
import matplotlib
matplotlib.use('Agg')  # don't need display
import matplotlib.pyplot as plt
from astropy.visualization import astropy_mpl_style
plt.style.use(astropy_mpl_style)
# python timezones
import pytz
import glob
import shutil

# import classes from the chultun module
from chultun import Target  # name, ra, dec
from chultun import Stack  # exposure, filter, binning, count
from chultun import Sequence  # stacks, repeat
from chultun import Observatory  # code, latitude, longitude, altitude
from chultun import Observation  # target, sequence
from chultun import Scheduler  # observatory, observations
from chultun import Telescope  # the telescope commands


def runSubprocess(command_array, simulate=False, communicate=True, timeout=0):
    # command array is array with command and all required parameters
    if simulate:
        logme('Simulating subprocess "%s".' % (command_array))
        return ('', 0, 0)
    try:
        if timeout > 0:
            command_array = ['timeout', str(timeout)] + command_array
        sp = subprocess.Popen(
            command_array, stderr=subprocess.PIPE, stdout=subprocess.PIPE)
        logme('Running subprocess ("%s" %s)...' %
              (' '.join(command_array), sp.pid))
        sp.wait()
        if communicate:
            output, error = sp.communicate(b'\n\n')
            if error:
                logme('Error. Process (%s) reported error.' %
                      (command_array))
            return (output, error, sp.pid)
        else:  # some processes, like keepopen, hang forever with .communicate()
            return ('', '', 0)
    except Exception as e:
        logme(traceback.format_exception(*sys.exc_info()))
        return ('', 'Unknown error.', 0)


def getStats(command, user):
    logme('Retrieving telescope statistics...')

    send_message('Retrieving telescope statistics...')
    send_message(
        'This may take several minutes. Remember, Tikal was not built in a day!')

    (output, error, pid) = runSubprocess(
        ['/home/mcnowinski/anaconda2/bin/python', '/home/mcnowinski/giterdone/seo/itzamna/unlogger.py'], simulate)

    stats = open('/home/mcnowinski/giterdone/seo/itzamna/stats.txt', 'w')
    stats.write(output)
    stats.close()

    send_file('/home/mcnowinski/giterdone/seo/itzamna/stats.txt', 'Metric')
    send_file('/home/mcnowinski/giterdone/seo/itzamna/cloud.png', 'Clouds')
    send_file('/home/mcnowinski/giterdone/seo/itzamna/pointing.png', 'Pointing')


def getObservability(command, user):
    match = re.search('^\\\\(plot)\\s?([0-9]+)?', command)
    if(match):
        if match.group(2):
            object_index = int(match.group(2).strip())
        else:
            logme('No object index provided in \\plot command. Choosing first object...')
            object_index = 1
    else:
        logme('Error. Unexpected command format (%s).' % command)
        return

    # if we don't have any objects from \find, go home
    if len(objects) <= 0:
        send_message(
            'Itzamna does not remember any objects. Please run `\\find <object>` to remind me.')
        return

    # if specified index does not match our object list, go home
    if object_index <= 0 or object_index > len(objects):
        send_message('Itzamna only remembers %d object(s):' % len(objects))
        report = ''
        index = 1
        # calculate local time of observatory
        # + object_observer_utc_offset
        object_observer_now = Time(datetime.datetime.utcnow(), scale='utc')
        for object in objects:
            # create SkyCoord instance from RA and DEC
            c = SkyCoord(object['RA'], object['DEC'], unit="deg")
            # transform RA,DEC to alt, az for this object from the observatory
            altaz = c.transform_to(
                AltAz(obstime=object_observer_now, location=object_observer))
            report += '%d.\t%s object (%s) found at RA=%s, DEC=%s, ALT=%f, AZ=%f.\n' % (
                index, object['type'], object['name'], object['RA'], object['DEC'], altaz.alt.degree, altaz.az.degree)
            index += 1
        send_message(report)
        return

    # select specific object from the find list
    object = objects[object_index-1]

    send_message(
        'Itzamna is calculating when "%s" is observable from your location. Please wait...' % object['name'])

    logme('Checking observability of %s (%d)...' %
          (object['name'], object_index))

    # calculate local time of observatory
    #observatory_utc_offset = int(datetime.datetime.now(pytz.timezone(observatory_tz)).strftime('%z'))/100.0*u.hour
    observer_now = Time(datetime.datetime.utcnow(),
                        scale='utc')  # + observatory_utc_offset
    observer_now_sidereal = observer_now.sidereal_time(
        'mean', observatory_lon*u.deg)

    # for celestial objects, this is easy peasy
    if object['type'] == 'Celestial':
        c = SkyCoord(object['RA'], object['DEC'], unit="deg")
        # look 24 hours into the future
        delta_now = numpy.linspace(0, 24, 1000)*u.hour
        times_now_to_tomorrow = observer_now + delta_now
        frame_now_to_tomorrow = AltAz(
            obstime=times_now_to_tomorrow, location=object_observer)
        object_altaz_now_to_tomorrow = c.transform_to(frame_now_to_tomorrow)
        sun_altaz_now_to_tomorrow = get_sun(
            times_now_to_tomorrow).transform_to(frame_now_to_tomorrow)
        #plt.plot(delta_now, sun_altaz_now_to_tomorrow.alt, color='r', label='Sun')
        delta_now_ha = numpy.linspace(0, 24, 1000)*u.hourangle
        times_now_to_tomorrow_sidereal = observer_now_sidereal + delta_now_ha
        #print times_now_to_tomorrow_sidereal
        hour_angle = numpy.empty((0))
        for ha in (times_now_to_tomorrow_sidereal-c.ra).hour:
            if ha >= 18:
                hour_angle = numpy.append(hour_angle, ha-24)
            elif ha <= -18:
                hour_angle = numpy.append(hour_angle, ha+24)
            else:
                hour_angle = numpy.append(hour_angle, ha)
        plt.scatter(delta_now, object_altaz_now_to_tomorrow.alt,
                    c=object_altaz_now_to_tomorrow.az, label=object['name'], lw=0, s=8,
                    cmap='viridis')
        plt.fill_between(delta_now.to('hr').value, 0, 90,
                         sun_altaz_now_to_tomorrow.alt < -0*u.deg, color='0.5', zorder=0)
        plt.fill_between(delta_now.to('hr').value, 0, 90,
                         sun_altaz_now_to_tomorrow.alt < -18*u.deg, color='k', zorder=0)
        plt.fill_between(delta_now.to('hr').value, 0, 90,
                         abs(hour_angle) <= 5.3, color='LightBlue', alpha=0.40, zorder=0)
        # plt.fill_between(delta_now.to('hr').value, 0, 90,
        #                 abs(hour_angle) <= 5.3, facecolor='None', hatch='\\', edgecolor='b', color='b')
        plt.colorbar().set_label('Azimuth [deg]')
        plt.legend(loc='upper left')
        plt.xlim(0, 24)
        plt.xticks(numpy.arange(13)*2)
        plt.ylim(0, 90)
        plt.xlabel('Hours [from now]')
        plt.ylabel('Altitude [deg]')
        # plt.show()
        plt.savefig('plot.png', bbox_inches='tight')
        plt.close()
        send_file('plot.png', 'Target (%s) Visibility' % object['name'])
        ##plt.plot(delta_now, sun_altaz_now_to_tomorrow.alt, color='r', label='Sun')
        # delta_now_ha = numpy.linspace(0, 24, 1000)*u.hourangle
        # times_now_to_tomorrow_sidereal = observer_now_sidereal + delta_now_ha
        ##print times_now_to_tomorrow_sidereal
        # hour_angle = numpy.empty((0))
        # for ha in (times_now_to_tomorrow_sidereal-c.ra).hour:
        # if ha >= 18:
        # hour_angle = numpy.append(hour_angle, ha-24)
        # elif ha <= -18:
        # hour_angle = numpy.append(hour_angle, ha+24)
        # else:
        # hour_angle = numpy.append(hour_angle, ha)
        # plt.scatter(delta_now, hour_angle,
        # c=object_altaz_now_to_tomorrow.az, label=object['name'], lw=0, s=8,
        # cmap='viridis')
        # plt.fill_between(delta_now.to('hr').value, -24, 24,
        # sun_altaz_now_to_tomorrow.alt < -0*u.deg, color='0.5', zorder=0)
        # plt.fill_between(delta_now.to('hr').value, -24, 24,
        # sun_altaz_now_to_tomorrow.alt < -18*u.deg, color='k', zorder=0)

        # plt.colorbar().set_label('Azimuth [deg]')
        # plt.legend(loc='upper left')
        # plt.xlim(0, 24)
        # plt.xticks(numpy.arange(13)*2)
        # plt.ylim(-24, 24)
        # plt.xlabel('Hours [from now]')
        # plt.ylabel('Hour Angle [hours]')
        # plt.show()
        # plt.savefig('plot.png', bbox_inches='tight')
        # plt.close()
        # send_file('plot.png', 'Target (%s) Visibility'%object['name'])
    elif object['type'] == 'Solar System':
        #send_message('Itzamna does not recognize the object type (%s).'%object['type'])
        # set time period to next 24 hours
        #observer_now = Time(datetime.datetime.utcnow(), scale='utc')
        start = observer_now.tt.datetime
        end = start+datetime.timedelta(days=1)
        # get object ephemerides for time range
        result = ch.query(object['id'], smallbody=False)
        result.set_epochrange(start.strftime(
            "%Y/%m/%d %H:%M"), end.strftime("%Y/%m/%d %H:%M"), '1m')
        result.get_ephemerides(observatory_code)
        # plot it
        delta_now = numpy.linspace(0, 24, 1441)*u.hour
        altitude = numpy.array(result['EL'])
        azimuth = numpy.array(result['AZ'])
        times_now_to_tomorrow = observer_now + delta_now
        #print times_now_to_tomorrow
        frame_now_to_tomorrow = AltAz(
            obstime=times_now_to_tomorrow, location=object_observer)
        sun_altaz_now_to_tomorrow = get_sun(
            times_now_to_tomorrow).transform_to(frame_now_to_tomorrow)
        plt.scatter(delta_now, altitude,
                    c=azimuth, label=object['name'], lw=0, s=8,
                    cmap='viridis')
        plt.fill_between(delta_now.to('hr').value, 0, 90,
                         sun_altaz_now_to_tomorrow.alt < -0*u.deg, color='0.5', zorder=0)
        plt.fill_between(delta_now.to('hr').value, 0, 90,
                         sun_altaz_now_to_tomorrow.alt < -18*u.deg, color='k', zorder=0)
        plt.colorbar().set_label('Azimuth [deg]')
        plt.legend(loc='upper left')
        plt.xlim(0, 24)
        plt.xticks(numpy.arange(13)*2)
        plt.ylim(0, 90)
        plt.xlabel('Hours [from now]')
        plt.ylabel('Altitude [deg]')
        plt.savefig('plot.png', bbox_inches='tight')
        plt.close()
        send_file('plot.png', 'Target (%s) Visibility' % object['name'])
    elif object['type'] == 'Satellite':
        #print object['id'], object['tle_line1'], object['tle_line2']
        sat_ephem = ephem.readtle(
            object['id'], object['tle_line1'], object['tle_line2'])
        # look 24 hours into the future
        delta_now = numpy.linspace(0, 24, 10000)*u.hour
        times_now_to_tomorrow = observer_now + delta_now
        frame_now_to_tomorrow = AltAz(
            obstime=times_now_to_tomorrow, location=object_observer)
        sun_altaz_now_to_tomorrow = get_sun(
            times_now_to_tomorrow).transform_to(frame_now_to_tomorrow)
        alt = []
        az = []
        for time_now_to_tomorrow in times_now_to_tomorrow:
            dt = time_now_to_tomorrow.tt.datetime
            sat_observer.date = dt
            sat_ephem.compute(sat_observer)
            #sat_coords = SkyCoord(ra='%s'%sat_ephem.ra, dec='%s'%sat_ephem.dec, unit=(u.hour, u.deg))
            alt.append(math.degrees(float(repr(sat_ephem.alt))))
            az.append(math.degrees(float(repr(sat_ephem.az))))
            #lat = math.degrees(float(repr(sat_ephem.sublat)))
            #lon = math.degrees(float(repr(sat_ephem.sublong)))
        #print numpy.mean(alt)
        plt.scatter(delta_now, alt,
                    c=az, label=object['name'], lw=0, s=8,
                    cmap='viridis')
        plt.fill_between(delta_now.to('hr').value, 0, 90,
                         sun_altaz_now_to_tomorrow.alt < -0*u.deg, color='0.5', zorder=0)
        plt.fill_between(delta_now.to('hr').value, 0, 90,
                         sun_altaz_now_to_tomorrow.alt < -18*u.deg, color='k', zorder=0)
        plt.colorbar().set_label('Azimuth [deg]')
        plt.legend(loc='upper left')
        plt.xlim(0, 24)
        plt.xticks(numpy.arange(13)*2)
        plt.ylim(0, 90)
        plt.xlabel('Hours [from now]')
        plt.ylabel('Altitude [deg]')
        # plt.show()
        plt.savefig('plot.png', bbox_inches='tight')
        plt.close()
        send_file('plot.png', 'Target (%s) Visibility' % object['name'])
    else:
        send_message(
            'Itzamna does not recognize the object type (%s).' % object['type'])


def doCrack(command, user):
    # this command requires that user has telescope locked
    if not lockedByYou(user):
        send_message('Please lock the telescope before calling this command.')
        return

    send_message("Itzamna is opening ('cracking') observatory. Please wait...")

    # reset the target name
    global target_name
    target_name = "unknown"

    (output, error, pid) = runSubprocess(['tin', 'interrupt'], simulate)
    logme(output)

    (output, error, pid) = runSubprocess(
        ['openup_nolock', 'nocloud'], simulate)
    logme(output)

    (output, error, pid) = runSubprocess(
        ['keepopen', 'maxtime=36000', 'slit'], simulate, False)
    logme(output)

    send_message("The observatory was successfully opened ('cracked')!")


def doSqueeze(command, user):
    # this command requires that user has telescope locked
    if not lockedByYou(user):
        send_message('Please lock the telescope before calling this command.')
        return

    send_message(
        "Itzamna is closing ('squeezing') the observatory. Please wait...")

    # reset the target name
    global target_name
    target_name = "unknown"

    (output, error, pid) = runSubprocess(['closedown_nolock'], simulate)
    if re.search('ERROR', output):
        send_message("Error! The observatory could not be closed!")
        # at least try to get slit closed!
        (output, error, pid) = runSubprocess(['tx', 'slit', 'close'], simulate)
        return
    logme(output)

    (output, error, pid) = runSubprocess(['tx', 'lock', 'clear'], simulate)
    if not re.search('done lock', output):
        logme('Error. Could not clear lock')

    (output, error, pid) = runSubprocess(['tin', 'resume'], simulate)
    logme(output)

    send_message("The observatory was successfully closed ('squoze')!")


def toStars(command, user):
    # are there any .fits images to send?
    fits = glob.glob(image_path+'*.fits')
    if(len(fits) <= 0):
        send_message('Itzamna does not have any recent images to send!')
        return
    else:
        send_message('Uploading %d image(s) to <%s|stars>. Itzamna is ready for your next command.' % (
            len(fits), stars_url))
    # we are going to put these in a folder corresponding to the datetime this command was run!
    # a bit different from how this usually works...

    dest_path = '%s%s/%s/%s' % (stars_image_path, datetime.datetime.utcnow().strftime(
        '%Y'), datetime.datetime.utcnow().strftime('%Y-%m-%d'), 'itzamna')
    # send_message(dest_path)
    # off they go!
    (output, error, pid) = runSubprocess(
        ['tostars', '%s*.fits' % image_path, '%s' % dest_path], simulate)
    if error == '':
        send_message(
            'Successfully uploaded %d image(s) to <%s|stars>!' % (len(fits), stars_url))
        # move images to archive
        files = glob.iglob(image_path+'*.fits')
        for file in files:
            if os.path.isfile(file):
                shutil.move(file, image_archive_path)
    else:
        send_message('Error. Image uploaded failed!')


def doImage(command, user):
    # this command requires that user has telescope locked
    if not lockedByYou(user):
        send_message('Please lock the telescope before calling this command.')
        return

    match = re.search(
        '^\\\\(image) ([0-9\\.]+) (0|1|2|3|4) (oiii|g\\-band|r\\-band|i\\-band|sii|clear|h\\-alpha)', command, re.IGNORECASE)
    if(match):
        exposure = match.group(2)
        binning = match.group(3)
        filter = match.group(4)
        send_message('Taking image (exposure=%s s, bin=%s, filter=%s). Please wait...' % (
            exposure, binning, filter))
    else:
        logme('Error. Unexpected command format (%s).' % command)
        return

    # IMAGE_FILENAME=${NAME}_${filter}_${EXPOSURE_SEC}s_bin${BINNING}_`date -u +"%y%m%d_%H%M%S"`__seo_${USER}_`printf "%04d" $COUNT`_RAW.fits
    fits = image_path + '%s_%s_%ss_bin%s_%s_%s_seo_%d_RAW.fits' % (
        target_name, filter, exposure, binning, datetime.datetime.utcnow().strftime('%y%m%d_%H%M%S'), 'itzamna', 0)
    fits = fits.replace(' ', '_')
    slackdebug('Taking image (%s)...' % (fits))
    #(output, error, pid) = runSubprocess(['pfilter', '%s' % filter], simulate)
    telescope.setFilter(filter)
    (output, error, pid) = runSubprocess(
        ['image', 'time=%s' % exposure, 'bin=%s' % binning, 'outfile=%s' % fits], simulate)
    if not error:
        send_message('Got image (%s).' % fits)
        slackpreview(fits)
    else:
        send_message('Error. Image command failed (%s).' % fits)

    (output, error, pid) = runSubprocess(['tx', 'track', 'on'], simulate)
    # done track ha=15.0410 dec=0.0000
    if not re.search('done track ha\\=[0-9\\+\\-\\.]+\\sdec\\=[0-9\\+\\-\\.]+', output):
        send_message('Error. Could not turn telescope tracking ON.')

# send alert message to slack


def doBias(command, user):
    # this command requires that user has telescope locked
    if not lockedByYou(user):
        send_message('Please lock the telescope before calling this command.')
        return

    match = re.search(
        '^\\\\(bias) (0|1|2|3|4)', command, re.IGNORECASE)
    if(match):
        exposure = '0.1'
        binning = match.group(2)
        filter = 'clear'
        send_message('Taking bias frame (bin=%s). Please wait...' % binning)
    else:
        logme('Error. Unexpected command format (%s).' % command)
        return
    # IMAGE_FILENAME=${NAME}_${filter}_${EXPOSURE_SEC}s_bin${BINNING}_`date -u +"%y%m%d_%H%M%S"`__seo_${USER}_`printf "%04d" $COUNT`_RAW.fits
    fits = image_path + '%s_%s_%ss_bin%s_%s_%s_seo_%d_RAW.fits' % (
        'bias', filter, exposure, binning, datetime.datetime.utcnow().strftime('%y%m%d_%H%M%S'), 'itzamna', 0)
    fits = fits.replace(' ', '_')
    slackdebug('Taking image (%s)...' % (fits))
    (output, error, pid) = runSubprocess(['pfilter', '%s' % filter], simulate)
    (output, error, pid) = runSubprocess(
        ['image', 'dark', 'time=%s' % exposure, 'bin=%s' % binning, 'outfile=%s' % fits], simulate)
    if not error:
        send_message('Got image (%s).' % fits)
        slackpreview(fits)
    else:
        send_message('Error. Image command failed (%s).' % fits)


def doDark(command, user):
    # this command requires that user has telescope locked
    if not lockedByYou(user):
        send_message('Please lock the telescope before calling this command.')
        return

    match = re.search(
        '^\\\\(dark) ([0-9\\.]+) (0|1|2|3|4)', command, re.IGNORECASE)
    if(match):
        exposure = match.group(2)
        binning = match.group(3)
        filter = 'h-alpha'
        send_message('Taking dark frame (exposure=%s, bin=%s). Please wait...' % (
            exposure, binning))
    else:
        logme('Error. Unexpected command format (%s).' % command)
        return

    # IMAGE_FILENAME=${NAME}_${filter}_${EXPOSURE_SEC}s_bin${BINNING}_`date -u +"%y%m%d_%H%M%S"`__seo_${USER}_`printf "%04d" $COUNT`_RAW.fits
    fits = image_path + '%s_%s_%ss_bin%s_%s_%s_seo_%d_RAW.fits' % (
        'dark', filter, exposure, binning, datetime.datetime.utcnow().strftime('%y%m%d_%H%M%S'), 'itzamna', 0)
    fits = fits.replace(' ', '_')
    slackdebug('Taking image (%s)...' % (fits))
    (output, error, pid) = runSubprocess(['pfilter', '%s' % filter], simulate)
    (output, error, pid) = runSubprocess(
        ['image', 'dark', 'time=%s' % exposure, 'bin=%s' % binning, 'outfile=%s' % fits], simulate)
    if not error:
        send_message('Got image (%s).' % fits)
        slackpreview(fits)
    else:
        send_message('Error. Image command failed (%s).' % fits)


# send alert message to slack
def slackdebugalert(msg):
    msg = datetime.datetime.utcnow().strftime('%m-%d-%Y %H:%M:%S ') + msg
    (output, error, pid) = runSubprocess(['slackalert', msg], simulate)

# send debug message to slack


def slackdebug(msg):
    msg = datetime.datetime.utcnow().strftime('%m-%d-%Y %H:%M:%S ') + msg
    (output, error, pid) = runSubprocess(['slackdebug', msg], simulate)

# send preview of fits image to Slack


def slackpreview(fits):
    (output, error, pid) = runSubprocess(
        ['stiffy', fits, 'image.tif'], simulate)
    (output, error, pid) = runSubprocess(
        ['convert', '-resize', '50%', '-normalize', '-quality', '75', 'image.tif', 'image.jpg'], simulate)
    (output, error, pid) = runSubprocess(
        ['slackpreview_itzamna', 'image.jpg', fits], simulate)
    (output, error, pid) = runSubprocess(
        ['rm', 'image.jpg', 'image.tif'], simulate)


def doShare(command, user):
    global share

    # this command requires that user has telescope locked
    # is telescope locked by us?
    (username, email) = lockedBy()
    if username != user['name']:
        send_message('Please lock the telescope before calling this command.')
        return

    match = re.search('^\\\\(share) (on|off)', command, re.IGNORECASE)
    if(match):
        toggle = match.group(2).lower()
        if(toggle == 'on'):
            share = True
            send_message('Your lock has been shared!')
        else:
            share = False
            send_message('Your lock has been un-shared!')
    else:
        logme('Error. Unexpected command format (%s).' % command)
        return


def doTrack(command, user):
    # this command requires that user has telescope locked
    if not lockedByYou(user):
        send_message('Please lock the telescope before calling this command.')
        return

    match = re.search('^\\\\(track) (on|off)', command, re.IGNORECASE)
    if(match):
        toggle = match.group(2).lower()
        (output, error, pid) = runSubprocess(['tx', 'track', toggle], simulate)
        # send_message(output)
        # done track ha=15.0410 dec=0.0000
        if not re.search('done track ha\\=[0-9\\+\\-\\.]+\\sdec\\=[0-9\\+\\-\\.]+', output):
            send_message(
                'Error. Could not set telescope tracking to %s.' % toggle)
        else:
            send_message('Telescope tracking is %s.' % toggle)
    else:
        logme('Error. Unexpected command format (%s).' % command)
        return


def doFocus(command, user):
    # this command requires that user has telescope locked
    if not lockedByYou(user):
        send_message('Please lock the telescope before calling this command.')
        return

    match = re.search('^\\\\(focus)\\s([0-9]+)', command, re.IGNORECASE)
    if(match):
        position = match.group(2)
        (output, error, pid) = runSubprocess(
            ['tx', 'focus', 'pos=%s' % position], simulate)
        # send_message(output)
        # done focus pos=4854
        if not re.search('done focus pos\\=([0-9]+)', output):
            send_message('Error. Could not set focus to %s.' % position)
        else:
            send_message('Focus position is %s.' % match.group(2))
    else:
        logme('Error. Unexpected command format (%s).' % command)
        return


def doOffset(command, user):
    # this command requires that user has telescope locked
    if not lockedByYou(user):
        send_message('Please lock the telescope before calling this command.')
        return

    match = re.search(
        '^\\\\(nudge)\\s([\\-\\.0-9]+)\\s([\\-\\.0-9]+)', command, re.IGNORECASE)
    if(match):
        offset_ra = float(match.group(2))  # in arcmin (of hour)
        offset_dec = float(match.group(3))  # in arcmin (of degree)
        if offset_ra > 60 or offset_ra < -60 or offset_dec > 60 or offset_dec < -60:
            send_message(
                'Error. Valid values of RA and DEC offsets are +/- 0 to 60 arcmin.')
            return
        # convert to degrees
        offset_ra_deg = (offset_ra/60.0)*(360.0/24.0)
        offset_dec_deg = offset_dec/60.0
        (output, error, pid) = runSubprocess(
            ['tx', 'offset', 'ra=%f' % offset_ra_deg, 'dec=%f' % offset_dec_deg], simulate)
        if not re.search('done offset', output):
            send_message("Error. Could not set offset to dRA=%.2f', dDEC=%.2f'." % (
                offset_ra, offset_dec))
        else:
            send_message("Telescope pointing successfully offset (dRA=%.2f', dDEC=%.2f')." % (
                offset_ra, offset_dec))
    else:
        logme('Error. Unexpected command format (%s).' % command)
        return


def doHomer(command, user):
    # this command requires that user has telescope locked
    if not lockedByYou(user):
        send_message('Please lock the telescope before calling this command.')
        return

    # close observatory before doing tuneup
    doSqueeze('', user)

    send_message('Telescope control system calibration started.')

    send_message('Homing HA drive...')
    (output, error, pid) = runSubprocess(['tx', 'home', 'ha'], simulate)

    send_message('Homing DEC drive...')
    (output, error, pid) = runSubprocess(['tx', 'home', 'dec'], simulate)

    send_message('Homing DOME drive...')
    (output, error, pid) = runSubprocess(['tx', 'home', 'domer'], simulate)
    (output, error, pid) = runSubprocess(['tx', 'home', 'domel'], simulate)

    send_message('Telescope control system calibration complete.')


def doPinpointByObjectNum(command, user):
    match = re.search('^\\\\(pinpoint)\\s?([0-9]+)?', command)
    if(match):
        if match.group(2):
            object_index = int(match.group(2).strip())
        else:
            logme(
                'No object index provided in \\pinpoint command. Choosing first object...')
            object_index = 1
    else:
        logme('Error. Unexpected command format (%s).' % command)
        return

    # this command requires that user has telescope locked
    if not lockedByYou(user):
        send_message('Please lock the telescope before calling this command.')
        return

    # slit should be open
    if getSlit('', user) != 'open':
        send_message(
            'Slit is not open. Please run \\crack to open observatory.')
        return

    # if we don't have any objects from \find, go home
    if len(objects) <= 0:
        send_message(
            'Itzamna does not remember any objects. Please run `\\find <object>` to remind me.')
        return

    # if specified index does not match our object list, go home
    if object_index <= 0 or object_index > len(objects):
        send_message('Itzamna only remembers %d object(s):' % len(objects))
        report = ''
        index = 1
        # calculate local time of observatory
        # + object_observer_utc_offset
        object_observer_now = Time(datetime.datetime.utcnow(), scale='utc')
        for object in objects:
            # create SkyCoord instance from RA and DEC
            c = SkyCoord(object['RA'], object['DEC'], unit="deg")
            # transform RA,DEC to alt, az for this object from the observatory
            altaz = c.transform_to(
                AltAz(obstime=object_observer_now, location=object_observer))
            ra = Angle('%fd' % object['RA']).to_string(unit=u.hour, sep=':')
            dec = Angle('%fd' % object['DEC']).to_string(
                unit=u.degree, sep=':')
            report += '%d.\t%s object (%s) found at RA=%s, DEC=%s, ALT=%f, AZ=%f.\n' % (
                index, object['type'], object['name'], ra, dec, altaz.alt.degree, altaz.az.degree)
            index += 1
        send_message(report)
        return

    # select specific object from the find list
    object = objects[object_index-1]

    send_message(
        'Itzamna is pinpointing the telescope to "%s". Please wait...' % object['name'])
    logme('Pinpointing the telescope to %s (%d)...' %
          (object['name'], object_index))

    # reset the target name
    global target_name
    target_name = re.sub('[^A-Za-z0-9]', '_', object['name'])

    (output, error, pid) = runSubprocess(['tx', 'track', 'on'], simulate)
    # send_message(output)
    if not re.search('done track ha\\=[0-9\\+\\-\\.]+\\sdec\\=[0-9\\+\\-\\.]+', output):
        send_message(
            'Error. Could not enable telescope tracking (%s).' % output)
        return

    ra = Angle('%fd' % object['RA']).to_string(unit=u.hour, sep=':')
    dec = Angle('%fd' % object['DEC']).to_string(unit=u.degree, sep=':')
    (output, error, pid) = runSubprocess(
        ['pinpoint', '%s' % ra, '%s' % dec], simulate, True, 300)
    # done point move=62.455 dist=0.0031
    # send_message(output)
    if not re.search('BAM', output):
        send_message('Error. Could not pinpoint the telescope (%s).' % output)
        return

    send_message('Telescope successfully pinpointed to "%s" (RA=%s, DEC=%s).' % (
        object['name'], ra, dec))


def doPinpointByRaDec(command, user):
    match = re.search(
        '^\\\\(pinpoint) ([0-9\\:\\-\\+\\.]+) ([0-9\\:\\-\\+\\.]+)', command)
    if(match):
        ra = match.group(2)
        dec = match.group(3)
    else:
        logme('Error. Unexpected command format (%s).' % command)
        return

    # this command requires that user has telescope locked
    if not lockedByYou(user):
        send_message('Please lock the telescope before calling this command.')
        return

    # slit should be open
    if getSlit('', user) != 'open':
        send_message(
            'Slit is not open. Please run \\crack to open observatory.')
        return

    # ensure these are valid coordinates
    try:
        ra_decimal = Angle(ra + '  hours')
        dec_decimal = Angle(dec + '  degrees')
    except:
        send_message('Invalid RA/DEC coordinates.')
        return

    send_message('Itzamna is pinpointing the telescope. Please wait...')

    # regex to format RA/dec for filename
    ra = re.sub('^(\d{1,2}):(\d{2}):(\d{2}).+', r'\1h\2m\3s', ra)
    dec = re.sub('(\d{1,2}):(\d{2}):(\d{2}).+', r'\1d\2m\3s', dec)

    # reset the target name
    global target_name
    target_name = "%s%s" % (ra, dec)

    (output, error, pid) = runSubprocess(['tx', 'track', 'on'], simulate)
    # send_message(output)
    if not re.search('done track ha\\=[0-9\\+\\-\\.]+\\sdec\\=[0-9\\+\\-\\.]+', output):
        send_message(
            'Error. Could not enable telescope tracking (%s).' % output)
        return

    (output, error, pid) = runSubprocess(
        ['pinpoint', '%s' % ra, '%s' % dec], simulate, True, 300)
    # done point move=62.455 dist=0.0031
    # send_message(output)
    if not re.search('BAM', output):
        send_message('Error. Could not pinpoint the telescope (%s).' % output)
        return

    send_message(
        'Telescope successfully pinpointed to RA=%s, DEC=%s.' % (ra, dec))


def doPointByObjectNum(command, user):
    match = re.search('^\\\\(point)\\s?([0-9]+)?', command)
    if(match):
        if match.group(2):
            object_index = int(match.group(2).strip())
        else:
            logme('No object index provided in \\point command. Choosing first object...')
            object_index = 1
    else:
        logme('Error. Unexpected command format (%s).' % command)
        return

    # this command requires that user has telescope locked
    if not lockedByYou(user):
        send_message('Please lock the telescope before calling this command.')
        return

    # if we don't have any objects from \find, go home
    if len(objects) <= 0:
        send_message(
            'Itzamna does not remember any objects. Please run `\\find <object>` to remind me.')
        return

    # if specified index does not match our object list, go home
    if object_index <= 0 or object_index > len(objects):
        send_message('Itzamna only remembers %d object(s):' % len(objects))
        report = ''
        index = 1
        # calculate local time of observatory
        # + object_observer_utc_offset
        object_observer_now = Time(datetime.datetime.utcnow(), scale='utc')
        for object in objects:
            # create SkyCoord instance from RA and DEC
            c = SkyCoord(object['RA'], object['DEC'], unit="deg")
            # transform RA,DEC to alt, az for this object from the observatory
            altaz = c.transform_to(
                AltAz(obstime=object_observer_now, location=object_observer))
            ra = Angle('%fd' % object['RA']).to_string(unit=u.hour, sep=':')
            dec = Angle('%fd' % object['DEC']).to_string(
                unit=u.degree, sep=':')
            report += '%d.\t%s object (%s) found at RA=%s, DEC=%s, ALT=%f, AZ=%f.\n' % (
                index, object['type'], object['name'], ra, dec, altaz.alt.degree, altaz.az.degree)
            index += 1
        send_message(report)
        return

    # select specific object from the find list
    object = objects[object_index-1]

    send_message(
        'Itzamna is pointing the telescope to "%s". Please wait...' % object['name'])
    logme('Pointing the telescope to %s (%d)...' %
          (object['name'], object_index))

    # reset the target name
    global target_name
    target_name = re.sub('[^A-Za-z0-9]', '_', object['name'])

    (output, error, pid) = runSubprocess(['tx', 'track', 'on'], simulate)
    if not re.search('done track ha\\=[0-9\\+\\-\\.]+\\sdec\\=[0-9\\+\\-\\.]+', output):
        send_message(
            'Error. Could not enable telescope tracking (%s).' % output)
        return

    ra = Angle('%fd' % object['RA']).to_string(unit=u.hour, sep=':')
    dec = Angle('%fd' % object['DEC']).to_string(unit=u.degree, sep=':')
    (output, error, pid) = runSubprocess(
        ['tx', 'point', 'ra=%s' % ra, 'dec=%s' % dec], simulate)
    # done point move=62.455 dist=0.0031
    if not re.search('done point move\\=[0-9\\+\\-\\.]+\\sdist\\=[0-9\\+\\-\\.]+', output):
        send_message('Error. Could not point the telescope (%s).' % output)
        return

    send_message('Telescope successfully pointed to "%s" (RA=%s, DEC=%s).' % (
        object['name'], ra, dec))


def doPointByRaDec(command, user):
    match = re.search(
        '^\\\\(point) ([0-9\\:\\-\\+\\.]+) ([0-9\\:\\-\\+\\.]+)', command)
    if(match):
        ra = match.group(2)
        dec = match.group(3)
    else:
        logme('Error. Unexpected command format (%s).' % command)
        return

    # this command requires that user has telescope locked
    if not lockedByYou(user):
        send_message('Please lock the telescope before calling this command.')
        return

    # ensure these are valid coordinates
    try:
        ra_decimal = Angle(ra + '  hours')
        dec_decimal = Angle(dec + '  degrees')
    except:
        send_message('Invalid RA/DEC coordinates.')
        return

    send_message('Itzamna is pointing the telescope. Please wait...')

    # regex to format RA/dec for filename
    ra = re.sub('^(\d{1,2}):(\d{2}):(\d{2}).+', r'\1h\2m\3s', ra)
    dec = re.sub('(\d{1,2}):(\d{2}):(\d{2}).+', r'\1d\2m\3s', dec)

    # reset the target name
    global target_name
    target_name = "%s%s" % (ra, dec)

    (output, error, pid) = runSubprocess(['tx', 'track', 'on'], simulate)
    if not re.search('done track ha\\=[0-9\\+\\-\\.]+\\sdec\\=[0-9\\+\\-\\.]+', output):
        send_message(
            'Error. Could not enable telescope tracking (%s).' % output)
        return

    (output, error, pid) = runSubprocess(
        ['tx', 'point', 'ra=%s' % ra, 'dec=%s' % dec], simulate)
    # returns done point move=62.455 dist=0.0031
    if not re.search('done point move\\=[0-9\\+\\-\\.]+\\sdist\\=[0-9\\+\\-\\.]+', output):
        send_message('Error. Could not point the telescope (%s).' % output)
        return

    send_message(
        'Telescope successfully pointed to RA=%s, DEC=%s.' % (ra, dec))


def getObject(command, user):
    global objects
    match = re.search('^\\\\(find) ([a-zA-Z0-9\\s\\-\\+\\*]+)', command)
    if(match):
        lookup = match.group(2)
    else:
        logme('Error. Unexpected command format (%s).' % command)
        return

    send_message(
        'Itzamna is searching the cosmos for "%s". Please wait...' % lookup)
    objects = []  # all objects found
    # max match results per type
    max_results = 50
    #
    # search celestial objects using SkyCoord
    #
    #found_celestial_object = False
    try:
        # add flux data to results
        Simbad.add_votable_fields('fluxdata(V)')
        result_table = Simbad.query_object(lookup.upper().replace('*', ''))
        if len(result_table) > max_results:
            send_message(
                'Exceeded maximum celestial matches. Will only return %d celestial results.' % (max_results))
        count = 0
        for row in range(0, len(result_table)):
            if count >= max_results:
                break
            count += 1
            objects.append({'type': 'Celestial', 'id': result_table['MAIN_ID'][row], 'name': result_table['MAIN_ID'][row].replace(' ', ''), 'RA': Angle(
                result_table['RA'][row].replace(' ', ':') + ' hours').degree, 'DEC': Angle(result_table['DEC'][row].replace(' ', ':') + ' degrees').degree, 'VMAG': result_table['FLUX_V'][row]})
    except:
        pass
    send_message('Found %d celestial match(es) for "%s".' %
                 (len(objects), lookup))
    #
    # search solar system small bodies using JPL HORIZONS
    #
    # list of matches
    object_names = []
    # set to * to make the searches wider by default
    suffix = ''
    # two passes, one for major (and maybe small) and one for (only) small bodies
    lookups = [lookup + suffix, lookup + suffix + ';']
    for repeat in range(0, 2):
        # user JPL Horizons batch to find matches
        f = urllib2.urlopen('https://ssd.jpl.nasa.gov/horizons_batch.cgi?batch=l&COMMAND="%s"' %
                            urllib2.quote(lookups[repeat].upper()))
        output = f.read()  # the whole enchilada
        lines = output.splitlines()  # line by line
        # no matches? go home
        if re.search('No matches found', output):
            logme('No matches found in JPL Horizons for %s.' %
                  lookups[repeat].upper())
        elif re.search('Target body name:', output):
            logme('Single match found in JPL Horizons for %s.' %
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
                    logme('Error. Could not parse id for single match major body (%s).' %
                          lookups[repeat].upper().replace(suffix, ''))
            else:
                # user search term is unique, so use it!
                object_names.append(
                    lookups[repeat].upper().replace(suffix, ''))
        elif repeat == 1 and re.search('Matching small-bodies', output):
            logme('Multiple small bodies found in JPL Horizons for %s.' %
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
                    # add semicolon for small body lookups
                    object_names.append(record_number + ';')
            # check our parse job
            match = re.search('(\\d+) matches\\.', output)
            if match:
                if int(match.group(1)) != match_count:
                    logme('Multiple JPL small body parsing error!')
                else:
                    logme('Multiple JPL small body parsing successful!')
        elif repeat == 0 and re.search('Multiple major-bodies', output):
            logme('Multiple major bodies found in JPL Horizons for %s.' %
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
                        # NO semicolon for major body lookups
                        object_names.append(record_number)
            # check our parse job
            match = re.search('Number of matches =([\\s\\d]+).', output)
            if match:
                if int(match.group(1)) != match_count:
                    logme('Multiple JPL major body parsing error!')
                else:
                    logme('Multiple JPL major body parsing successful!')
    start = datetime.datetime.utcnow()
    end = start+datetime.timedelta(seconds=60)
    send_message('Found %d solar system match(es) for "%s".' %
                 (len(object_names), lookup))
    if len(object_names) > max_results:
        send_message(
            'Exceeded maximum solar system matches. Will only return %d solar system results.' % (max_results))
    count = 0
    if len(object_names):
        send_message(
            'Calculating current sky positions of solar system match(es)...')
    for object_name in object_names:
        if count >= max_results:
            break
        count += 1
        try:
            result = ch.query(object_name.upper(), smallbody=False)
            result.set_epochrange(start.strftime(
                "%Y/%m/%d %H:%M"), end.strftime("%Y/%m/%d %H:%M"), '1m')
            result.get_ephemerides(observatory_code)
            objects.append({'type': 'Solar System', 'id': object_name.upper(
            ), 'name': result['targetname'][0], 'RA': result['RA'][0], 'DEC': result['DEC'][0], 'VMAG': result['V'][0]})
        except:
            pass
    #
    # search satellites
    #
    num_sat_matches = 0
    for sat in norad_sat_db:
        got_match = False
        if lookup.find('*') >= 0:
            got_match = (sat[0].find(lookup.upper().replace('*', '')) >= 0)
        else:
            got_match = (sat[0] == lookup.upper())
        if got_match:
            if num_sat_matches >= max_results:
                break
            num_sat_matches += 1
            sat_name = sat[0]
            sat_tle_line1 = sat[1]
            sat_tle_line2 = sat[2]
            sat_ephem = ephem.readtle(sat_name, sat_tle_line1, sat_tle_line2)
            sat_observer.date = datetime.datetime.utcnow()
            sat_ephem.compute(sat_observer)
            sat_coords = SkyCoord(ra='%s' % sat_ephem.ra, dec='%s' %
                                  sat_ephem.dec, unit=(u.hour, u.deg))
            objects.append({'type': 'Satellite', 'tle_line1': sat_tle_line1, 'tle_line2': sat_tle_line2, 'id': sat_name,
                            'name': sat_name, 'RA': math.degrees(float(repr(sat_ephem.ra))), 'VMAG': '', 'DEC': math.degrees(float(repr(sat_ephem.dec)))})

    send_message('Found %d satellite match(es) for "%s".' %
                 (num_sat_matches, lookup))
    if num_sat_matches > max_results:
        send_message(
            'Exceeded maximum satellite matches. Will only return %d satellite results.' % (max_results))
    #
    # process total search restults
    # calc current alt and az, then send information to Slack
    if len(objects) > 0:
        report = ''
        index = 1
        # calculate local time of observatory
        object_observer_now = Time(datetime.datetime.utcnow(), scale='utc')
        #print object_observer_now
        for object in objects:
            # create SkyCoord instance from RA and DEC
            c = SkyCoord(object['RA'], object['DEC'], unit="deg")
            # transform RA,DEC to alt, az for this object from the observatory
            altaz = c.transform_to(
                AltAz(obstime=object_observer_now, location=object_observer))
            ra = Angle('%fd' % object['RA']).to_string(unit=u.hour, sep=':')
            dec = Angle('%fd' % object['DEC']).to_string(
                unit=u.degree, sep=':')
            report += '%d.\t%s object (%s) found at RA=%s, DEC=%s, ALT=%f, AZ=%f, VMAG=%s.\n' % (
                index, object['type'], object['name'], ra, dec, altaz.alt.degree, altaz.az.degree, object['VMAG'])
            index += 1
        report = report.replace("--", "N/A")
        send_message(report)
    else:
        send_message(
            'Sorry, Itzamna knows all but *still* could not find "%s".' % lookup)


def lockedBy():
    logme('Checking to see if the telescope is locked...')

    locked_by = (None, None)
    (output, error, pid) = runSubprocess(['tx', 'lock'], simulate)
    # done lock user=mcnowinski email=mcnowinski@gmail.com phone=7032869140 comment=slac timestamp=2017-02-10T20:32:03Z
    match = re.search('^done lock user=(\S+) email=(\S+)', output)
    if(match):
        locked_by = (match.group(1), match.group(2))
        logme('Telescope is currently locked by %s (%s).' %
              (match.group(1), match.group(2)))
    return locked_by


def lockedByYou(user):
    if simulate:
        return True
    # check to make sure the telescope is locked by us!
    (username, email) = lockedBy()
    if username != user['name'] and share == False:
        return False
    else:
        return True


def doLock(command, user):
    logme('Locking the telescope...')

    # check to make sure the telescope isn't already locked
    (username, email) = lockedBy()
    if username:
        send_message('The telescope is already locked by %s (%s)!' %
                     (username, email))
        return False

    username = user['name']
    email = user['profile']['email']
    # only grab numbers from phone number
    num = re.compile(r'[^\d]+')
    try:
        phone = num.sub('', user['profile']['phone'])
    except:
        phone = '0000000000'
    comment = 'itzamna'
    timestamp = datetime.datetime.utcnow().strftime('%Y-%m-%dT%H:%M:%SZ')

    (output, error, pid) = runSubprocess(
        ['tx', 'lock', 'user=%s' % username, 'email=%s' % email, 'phone=%s' % phone, 'comment=%s' % comment], simulate)
    match = re.search('^done lock user=(\S+) email=(\S+)', output)
    if(match):
        send_message('Telescope successfully locked!')
    else:
        logme('Error. Telescope could not be locked (%s)' % output)
        send_message('Telescope could *not* be locked!')


def doUnLock(command, user):
    global share

    logme('Unlocking the telescope...')

    # check to make sure the telescope is locked and if so, locked by us!
    (username, email) = lockedBy()
    # is telescope locked?
    if not username:
        send_message("The telescope is not locked.")
        return False
    # is telescope locked by us?
    if username != user['name']:
        logme("Warning. Can't unlock the telescope. The telescope is already locked by %s (%s)!" % (
            username, email))
        send_message('The telescope is already locked by %s (%s)!' %
                     (username, email))
        return False

    # clear the lock
    (output, error, pid) = runSubprocess(['tx', 'lock', 'clear'], simulate)
    match = re.search('^done lock$', output)
    if(match):
        send_message('Telescope successfully unlocked!')
        share = False
    else:
        logme('Error. Telescope could not be unlocked (%s)' % output)
        send_message('Telescope could *not* be unlocked!')


def getSlit(command, user):
    logme('Retrieving dome slit status...')

    status = 'unknown'

    (output, error, pid) = runSubprocess(['tx', 'slit'], simulate)
    # done slit slit=closed
    match = re.search('^done slit slit=([a-zA-Z]+)', output)
    if(match):
        status = match.group(1)
    else:
        send_message(
            'Error. Command (%s) did not return a valid response.' % command)
        logme('Error. Command (%s) did not return a valid response (%s).' %
              (command, output))

    return status


def getSun(command, user):
    logme('Retrieving current sun information...')

    (output, error, pid) = runSubprocess(['sun'], simulate)
    # 21:30:51.07 -14:42:54.0 2017.107 sun alt=35.8
    match = re.search(
        '^([\\-\\+0-9\\:\\.]+) ([\\-\\+0-9\\:\\.]+) ([\\-\\+0-9\\.]+) sun alt=([\\-\\+0-9\\.]+)', output)
    if(match):
        send_message('Sun: RA=%s, DEC=%s, Alt=%s deg' %
                     (match.group(1), match.group(2), match.group(4)))
        send_message('\n')
    else:
        send_message(
            'Error. Command (%s) did not return a valid response.' % command)
        logme('Error. Command (%s) did not return a valid response (%s).' %
              (command, output))


def getClouds(command, user):
    logme('Retrieving current cloud cover...')

    (output, error, pid) = runSubprocess(['tx', 'taux'], simulate)
    # done taux ovolts=3.048 irvolts=0.199 cloud=0.26 rain=0 dew=2.97
    match = re.search('cloud=(\S+)', output)
    if(match):
        send_message('Cloud cover is %d%%.' % (int(float(match.group(1))*100)))
        send_message('\n')
    else:
        send_message(
            'Error. Command (%s) did not return a valid response.' % command)
        logme('Error. Command (%s) did not return a valid response (%s).' %
              (command, output))


def getFocus(command, user):
    logme('Retrieving current focus...')

    (output, error, pid) = runSubprocess(['tx', 'focus'], simulate)
    # done focus pos=4854
    match = re.search('pos=(\S+)', output)
    if(match):
        send_message('Focus position is %d.' % (int(match.group(1))))
        send_message('\n')
    else:
        send_message(
            'Error. Command (%s) did not return a valid response.' % command)
        logme('Error. Command (%s) did not return a valid response (%s).' %
              (command, output))


def getMoon(command, user):
    logme('Retrieving current moon information...')

    (output, error, pid) = runSubprocess(['moon'], simulate)
    # 07:38:11.68 +17:30:06.9 2017.107 moon alt=-21.1 phase=0.85 lunation=1164
    match = re.search(
        '^([\\-\\+0-9\\:\\.]+) ([\\-\\+0-9\\:\\.]+) ([\\-\\+0-9\\.]+) moon alt=([\\-\\+0-9\\.]+) phase=([\\-\\+0-9\\.]+) lunation=([\\-\\+0-9]+)', output)
    if(match):
        send_message('Moon: Phase=%d%%, RA=%s, DEC=%s, Alt=%s deg' % (
            int(float(match.group(5))*100.0), match.group(1), match.group(2), match.group(4)))
        send_message('\n')
    else:
        send_message(
            'Error. Command (%s) did not return a valid response.' % command)
        logme('Error. Command (%s) did not return a valid response (%s).' %
              (command, output))


def getWhere(command, user):
    logme('Retrieving the current telescope pointing information...')

    (output, error, pid) = runSubprocess(['tx', 'where'], simulate)
    # done where ra=05:25:25.11 dec=+38:17:17.0 equinox=2017.105 ha=0.010 secz=1.00 alt=90.0 az=265.1 slewing=0
    match = re.search(
        '^done where ra=([\\-\\+0-9\\:\\.]+) dec=([\\-\\+0-9\\:\\.]+) equinox=([\\-\\+0-9\\.]+) ha=([\\-\\+0-9\\.]+) secz=([\\-\\+0-9\\.]+) alt=([\\-\\+0-9\\.]+) az=([\\-\\+0-9\\.]+) slewing=([\\-\\+0-9\\.]+)', output)
    if(match):
        send_message('Telescope Pointing:')
        send_message('>RA: %s' % match.group(1))
        send_message('>DEC: %s' % match.group(2))
        send_message('>Alt: %s' % match.group(6))
        send_message('>Az: %s' % match.group(7))
        is_slewing = int(match.group(8))
        if is_slewing == 0:
            send_message('>Slewing? No')
        else:
            send_message('>Slewing? Yes')
        ra_decimal = Angle(match.group(1) + '  hours')
        dec_decimal = Angle(match.group(2) + '  degrees')
        url = 'http://server3.sky-map.org/imgcut?survey=DSS2&img_id=all&angle=0.5&ra=%f&de=%f&width=400&height=400&projection=tan&interpolation=bicubic&jpeg_quality=0.8&output_type=jpeg' % (
            ra_decimal.hour, dec_decimal.degree)
        send_message("", [{"image_url": "%s" %
                           url, "title": "Sky Position (DSS2):"}])
        send_message('\n')
    else:
        send_message(
            'Error. Command (%s) did not return a valid response.' % command)
        logme('Error. Command (%s) did not return a valid response (%s).' %
              (command, output))


def doWelcome():
    logme('Sending welcome message to Slack users...')

    send_message("", [{"image_url": "%s" % welcome_giphy_url,
                       "title": "Itzamna is here! Let your wishes be known..."}])
    # show help
    getHelp('\\help')

    lockedBy()

# get ClearDarkSky chart


def getClearDarkSky(command, user):
    logme('Retrieving the current Clear Sky charts for SEO...')

    dummy = ''.join(random.choice(string.ascii_uppercase + string.digits)
                    for _ in range(5))
    send_message("", [{"image_url": "http://www.cleardarksky.com/c/SonomaCAcsk.gif?dummy=%s" %
                       dummy, "title": "Lake Sonoma Clear Sky Chart"}])
    send_message("", [{"image_url": "http://www.cleardarksky.com/c/SmnCAcsk.gif?dummy=%s" %
                       dummy, "title": "Sonoma Clear Sky Chart"}])
    send_message("\n")


def getSkyCam(command, user):

    logme('Retrieving skycam image from SEO spacam...')

    (output, error, pid) = runSubprocess(['spacam'], False)
    send_file('spacam.jpg', 'SEO Spa-Cam in El Verano, CA')

    logme('Retrieving skycam images for sites near SEO...')

    dummy = ''.join(random.choice(string.ascii_uppercase + string.digits)
                    for _ in range(5))
    send_message("", [{"image_url": "http://icons.wunderground.com/webcamramdisk/c/v/cvogeo/1/current.jpg?%s" %
                       dummy, "title": "cvogeo's Webcam in Petaluma, CA"}])
    send_message("", [{"image_url": "http://icons.wunderground.com/webcamramdisk/w/u/WU_3508672/2/current.jpg?%s" %
                       dummy, "title": "WU_3508672's Webcam in Kenwood, CA"}])
    send_message("", [{"image_url": "http://icons.wunderground.com/webcamramdisk/l/p/lparkerwu66/1/current.jpg?%s" %
                       dummy, "title": "lparkerwu66's Webcam in Santa Rosa, CA"}])
    send_message("", [{"image_url": "http://icons.wunderground.com/webcamramdisk/j/w/JWPAGE/1/current.jpg?%s" %
                       dummy, "title": "JWPAGE's Webcam in Petaluma, CA"}])

    send_message("\n")
# get weather from Wunderground


def getForecast(command, user):
    logme('Retrieving the hourly forecast from wunderground.com...')

    #try:
    f = urllib2.urlopen('http://api.wunderground.com/api/%s/geolookup/hourly/q/pws:%s.json' %
                        (wunderground_token, wunderground_station))
    #except:
    #    send_message('Connection to wunderground failed.')
    #    return
    json_string = f.read()
    parsed_json = json.loads(json_string)
    hourly_forecasts = parsed_json['hourly_forecast']
    count = 0
    send_message("Weather Forecast:")
    for hourly_forecast in hourly_forecasts:
        count += 1
        if count > wunderground_max_forecast_hours:
            break
        send_message("", [{"image_url": "%s" % hourly_forecast['icon_url'], "title":"%s at %s:" % (
            hourly_forecast['condition'], hourly_forecast['FCTTIME']['pretty'])}])
    send_message("\n")
    f.close()

# get weather from Wunderground


def getWeather(command, user):
    logme('Retrieving the current weather conditions from wunderground.com...')

    # just in case wunderground is down...
    try:
        f = urllib2.urlopen('http://api.wunderground.com/api/%s/geolookup/conditions/q/pws:%s.json' %
                            (wunderground_token, wunderground_station))
    except:
        send_message("Weather forecast failed.")
        return
    json_string = f.read()
    parsed_json = json.loads(json_string)
    location = parsed_json['current_observation']['observation_location']['city']
    station = parsed_json['current_observation']['station_id']
    temp = parsed_json['current_observation']['temperature_string']
    weather = parsed_json['current_observation']['weather']
    rh = parsed_json['current_observation']['relative_humidity']
    wind = parsed_json['current_observation']['wind_string']
    wind_dir = parsed_json['current_observation']['wind_dir']
    wind_mph = parsed_json['current_observation']['wind_mph']
    dewpoint = parsed_json['current_observation']['dewpoint_string']
    icon_url = parsed_json['current_observation']['icon_url']
    last_update = parsed_json['current_observation']['observation_time']
    send_message("", [{"image_url": "%s" %
                       icon_url, "title": "Current Weather:"}])
    send_message(">%s" % (last_update))
    send_message(">Conditions: %s" % (weather))
    send_message(">Temperature: %s" % (temp))
    send_message(">Winds: %s" % (wind))
    send_message(">Humidity: %s" % (rh))
    send_message(">Local Station: %s (%s)" % (location, station))
    send_message("\n")
    f.close()


def getHelp(command, user=None):
    logme('Processing the "help" command...')

    # allow getHelp to be called by Itzamna
    user_name = 'Fear not, mortals'
    if user != None:
        user_name = user['profile']['display_name']
    send_message(user_name + ', here are some helpful tips:\n' +
                 '>Please report itzamna issues here: https://github.com/mcnowinski/seo/issues/new\n' +
                 '>A more detailed itzamna tutorial can be found here: https://stoneedgeobservatory.com/guide-to-using-itzamna/\n' +
                 '>`\\help` shows this message\n' +
                 '>`\\where` shows where the telescope is pointing\n' +
                 '>`\\weather` shows the current weather conditions\n' +
                 '>`\\forecast` shows the hourly weather forecast\n' +
                 '>`\\clouds` shows the current cloud cover\n' +
                 '>`\\focus` shows the current focus position\n' +
                 '>`\\focus <position>` sets the current focus position\n' + \
                 # '>`\\stats` shows the weekly telescope statistics\n' + \
                 '>`\\clearsky` shows the Clear Sky chart(s)\n' + \
                 '>`\\skycam` shows nearby skycam images\n' + \
                 '>`\\find <object>` finds <object> position in sky (add wildcard `*` to widen search)\n' + \
                 '>`\\plot <object#>` shows if/when <object> is observable (run `\\find` first!)\n' + \
                 '>`\\lock` locks the telescope\n' + \
                 '>`\\unlock` unlocks the telescope\n' + \
                 '>`\\crack` opens the observatory\n' + \
                 '>`\\squeeze` closes the observatory\n' + \
                 '>`\\share <on/off>` shares/unshares a locked telescope with others\n' + \
                 #'>`\\homer` re-homes the scope and dome (this will `\squeeze` the observatory!)\n'
                 '>`\\point <RA (hh:mm:ss.s)> <DEC (dd:mm:ss.s)>` or `\\point <object#>` points the telescope\n' + \
                 '>`\\pinpoint <RA (hh:mm:ss.s)> <DEC (dd:mm:ss.s)>` or `\\pinpoint <object#>` pinpoints the telescope\n' + \
                 # '>`\\track <on/off>` toggles telescope tracking\n' + \
                 # '>`\\nudge <dRA in arcmin> <dDEC in arcmin>` offsets the telescope pointing\n' + \
                 '>`\\image <exposure> <binning> <filter>` takes a picture\n' + \
                 '>`\\bias <binning>` takes a bias frame\n' + \
                 '>`\\dark <exposure> <binning>` takes a dark frame.\n' + \
                 '>`\\tostars` uploads recent images to <%s|stars> (run this command at the end of your session)\n' % stars_url
                 )
    send_message('\n')


def abort(msg):
    logme(msg)
    os.sys.exit(1)

# print and log messages


def logme(msg):
    # open log file
    log = open(log_fname, 'a+')
    dt = datetime.datetime.now().strftime("%Y/%m/%d %H:%M:%S:\t")
    log.write((dt + msg + "\n").encode('utf8'))
    log.close()
    print (dt + msg).encode('utf8')
    return

# send a message into the slack_channel


def send_message(msg, attachments=None):
    try:
        if slack_connected:
            sc.api_call(
                "chat.postMessage",
                channel=slack_channel,
                text=msg+'\n',
                username=bot_name,
                attachments=attachments
            )
            if len(msg.strip()) > 0:
                logme('Slack message: "%s"' % (msg.strip()))
                return
    except:
        pass
    logme('Error! Could not send message. Client is not connected.')


def send_file(path, title):
    if slack_connected:
        # curl -F file=@$file -F channels=$slackpreview_channel -F title="$title" -F token=$slackpreview_token $slackpreview_url > /dev/null 2>&1
        files = {'file': open(path, 'rb')}
        data = {'channels': slack_channel_name,
                'title': title, 'token': slack_token}
        url = slack_file_upload_url
        r = requests.post(url, files=files, data=data)
    else:
        logme('Error! Could not send file (%s). Client is not connected.' % path)

# get a list of slack users


def get_users():
    try:
        if slack_connected:
            result = sc.api_call("users.list")
            if 'members' in result:
                return result['members']
    except:
        pass
    logme('Error! Could not get user list. Client is not connected.')
    return []

# get a list of private and public channels


def get_channels():
    try:
        if slack_connected:
            channels = []
            # combined public and private channels
            result = sc.api_call("channels.list")
            if 'channels' in result:
                channels += result['channels']
            result = sc.api_call("groups.list")
            if 'groups' in result:
                channels += result['groups']
            if len(channels):
                return channels
    except:
        pass
    logme('Error! Could not get channel list. Client is not connected.')
    return []


def buildSatDatabase():
    global norad_sat_db
    for url in norad_sats_urls:
        # grab NORAD geosat data
        match = re.search('zip$', url)
        # if it's a zipped file, unzip first!
        if match:
            zipfile = ZipFile(StringIO(urllib2.urlopen(url).read()))
            sats = zipfile.open(zipfile.namelist()[0]).readlines()
        else:
            try:
                sats = urllib2.urlopen(url).readlines()
            except:
                logme('Error! Failed to open the satellite database.')
                return []
        # clean it up
        sats = [item.strip() for item in sats]
        # create an array of name, tle1, and tle2
        sats = [(str.upper(sats[i]), sats[i+1], sats[i+2])
                for i in xrange(0, len(sats)-2, 3)]
        # add sats to norad database
        norad_sat_db = numpy.concatenate((norad_sat_db, sats))
    logme('Read NORAD Two-Line Elements for %d satellite(s).' %
          len(norad_sat_db))


def ping():
    if slack_connected:
        try:
            data = json.dumps({"type": "ping"})
            sc.server.websocket.send(data)
            return True
        except:
            return False


def process_messages(msgs):
    global dt_last_message
    global dt_last_activity
    global slack_users
    for msg in msgs:
        dt_last_activity = datetime.datetime.now()
        # look for a msg from a user in the slack_channel, itzamna
        if 'type' in msg and 'channel' in msg and 'user' in msg:
            if msg['type'] == 'message' and msg['channel'] == slack_channel:
                dt = datetime.datetime.fromtimestamp(float(msg['ts']))
                # which user sent this message?
                # user_name='Unknown'
                user = None
                # update list of users, if get_users is successful (sometimes there are connection probs)
                latest_slack_users = get_users()
                if len(latest_slack_users) > 0:
                    slack_users = latest_slack_users
                for slack_user in slack_users:
                    if(slack_user['id'] == msg['user']):
                        user = slack_user
                        break
                if user == None:  # make a fake one
                    logme('User not found, using default values...')
                    user = {'name': 'Slack User', 'profile': {
                        'display_name': 'Slack User', 'email': 'user@slack.com', 'phone': '8675309'}}

                if(dt > dt_last_message):
                    dt_last_message = dt
                    text = msg['text'].strip()
                    match = re.search('^\\\\', text)
                    if match:
                        parse_command(text, user, dt)
                    else:
                        try:
                            logme('User %s sent text (%s) on %s.' % (
                                user['profile']['display_name'], msg['text'], dt_last_message.strftime("%Y/%m/%d @ %H:%M:%S")))
                        except:
                            pass
                    # is this a command, starts with \
                else:
                    logme('Warning! Ignoring old/duplicate message from #%s ("%s" from %s).' %
                          (slack_channel_name, msg['text'], user['profile']['display_name']))


def parse_command(text, user, dt):
    match = False
    for command in commands:
        # match the command first, then look for
        match = re.search(command[0], text, re.IGNORECASE)
        if match:  # compare with list of known commands, spelling and capitalization count!
            logme('%s sent command (%s).' %
                  (user['profile']['display_name'], text))
            # call associated function
            command[1](text, user)
            break
    if not match:  # did not recognize this command
        logme('%s sent unrecognized command (%s) on %s.' % (
            user['profile']['display_name'], text, dt.strftime("%Y/%m/%d @ %H:%M:%S")))
        send_message('%s, the almighty Itzamna does not recognize your command (%s).' % (
            user['profile']['display_name'], text))


def doSimulate():
    logme('Configuring simulate mode...')
    # use #dev for testing
    global slack_channel_name
    slack_channel_name = 'dev'


###############################
#CHANGE THESE VALUES AS NEEDED#
###############################
# run in simulate mode? restrict telescope commands
simulate = False
# log file
log_fname = 'itzamna.log'
# name of channel assigned to telescope interface
slack_channel_name = 'itzamna'
# how long to wait before successive slack reads
read_delay_s = 1
# how long to wait before trying to reconnect after connection fails or drops
reconnect_delay_s = 10
# name of this book
bot_name = 'Itzamna'
# specify a station close to SEO, e.g. LOLO Sonoma Farms
#wunderground_station = 'KCASONOM27'
#wunderground_station = 'KCASONOM51'
wunderground_station = 'KCASONOM64'
# how many hours of forecast should we show?
wunderground_max_forecast_hours = 12
# giphy shown when itzamna app is first started
welcome_giphy_url = 'http://www.nowinski.com/downloads/itzamna.gif'
# norad sat tle database urls
norad_sats_urls = [
    'http://www.celestrak.com/NORAD/elements/geo.txt',
    'http://www.celestrak.com/NORAD/elements/iridium.txt',
    'http://www.celestrak.com/NORAD/elements/iridium-NEXT.txt',
    'https://www.celestrak.com/NORAD/elements/tle-new.txt',
    'https://www.prismnet.com/~mmccants/tles/inttles.zip'
]
# make sure these folders exist!
image_path = '/home/mcnowinski/itzamna/images/'
image_archive_path = '/home/mcnowinski/itzamna/archive/'
# stars.uchicago.edu
stars_image_path = '/data/images/StoneEdge/0.5meter/'
###############################
#CHANGE THESE VALUES AS NEEDED#
###############################

# valid commands
# element 1 is a regular expression describing the command (and its parameters)
# element 2 is the function called when this command is received
commands = [

    ['^\\\\(help)', getHelp],
    #['^\\\\(weather)(\s)?(hourly)?', getWeather],
    ['^\\\\(weather)', getWeather],
    ['^\\\\(forecast)', getForecast],
    ['^\\\\(clearsky)', getClearDarkSky],
    ['^\\\\(skycam)', getSkyCam],
    ['^\\\\(where)', getWhere],
    ['^\\\\(sun)', getSun],
    ['^\\\\(moon)', getMoon],
    ['^\\\\(slit)', getSlit],
    #['^\\\\(test)', doTest],
    ['^\\\\(find) ([a-zA-Z0-9\\s\\-\\+\\*]+)', getObject],
    ['^\\\\(plot)\\s?([0-9]+)?', getObservability],
    ['^\\\\(stats)', getStats],
    ['^\\\\(point) ([0-9\\:\\-\\+\\.]+) ([0-9\\:\\-\\+\\.]+)', doPointByRaDec],
    ['^\\\\(point)\\s?([0-9]+)?', doPointByObjectNum],
    ['^\\\\(pinpoint) ([0-9\\:\\-\\+\\.]+) ([0-9\\:\\-\\+\\.]+)',
     doPinpointByRaDec],
    ['^\\\\(pinpoint)\\s?([0-9]+)?', doPinpointByObjectNum],
    ['^\\\\(track) (on|off)', doTrack],
    ['^\\\\(crack)', doCrack],
    ['^\\\\(squeeze)', doSqueeze],
    ['^\\\\(image) ([0-9\\.]+) (0|1|2|3|4) (oiii|g\\-band|r\\-band|i\\-band|sii|clear|h\\-alpha)', doImage],
    ['^\\\\(bias) (0|1|2|3|4)', doBias],
    ['^\\\\(dark) ([0-9\\.]+) (0|1|2|3|4)', doDark],
    ['^\\\\(lock)', doLock],
    ['^\\\\(share) (on|off)', doShare],
    ['^\\\\(unlock)', doUnLock],
    ['^\\\\(homer)', doHomer],
    ['^\\\\(clouds)', getClouds],
    ['^\\\\(focus)\\s([0-9]+)', doFocus],
    ['^\\\\(focus)', getFocus],
    ['^\\\\(nudge)\\s([\\-\\.0-9]+)\\s([\\-\\.0-9]+)', doOffset],
    ['^\\\\(tostars)', toStars]
]

# ensure slack token has been provided
if(len(sys.argv) < 3):
    abort('Error! Invalid command line arguments. Use "itzanma <SLACK_API_TOKEN> <WUNDERGROUND_API_TOKEN>".')

# the Slack api token
slack_token = sys.argv[1]
# the Slack url for file uploads
slack_file_upload_url = "https://slack.com/api/files.upload"
# the stars server URL (2018)
stars_url = 'http://stars.uchicago.edu/fitsview18/'
# the Wunderground api token
wunderground_token = sys.argv[2]
# track if slack client is connected
slack_connected = False
# list of slack users
slack_users = []
# list of slack channels
slack_channels = []
# id of slack channel assigned to telescope interface
slack_channel = None
# track time of last message received
# use to ignore old/duplcate messages
dt_last_message = dt_last_activity = datetime.datetime.now()
# observatory info
# seo
observatory_code = 'G52'
observatory_lat = 38.259  # deg
observatory_lon = -122.440  # deg
observatory_elev = 63.8  # m
observatory_tz = 'US/Pacific'  # for pytz
# satellite database; starts empty
norad_sat_db = numpy.array([]).reshape(0, 3)
# set pyephem observer
sat_observer = ephem.Observer()
sat_observer.lat, sat_observer.lon = '%f' % observatory_lat, '%f' % observatory_lon
# calc SkyCoord Earth location of observatory
object_observer = EarthLocation(
    lat=observatory_lat*u.deg, lon=observatory_lon*u.deg, height=observatory_elev*u.m)
# current object list
objects = []  # all objects found
# current target name
target_name = "unknown"
# lock sharing
share = False

if simulate:
    doSimulate()  # adjust configuration for simulated mode

logme('Starting itzamna Slack bot service...')

# build up a database of satellites from NORAD TLEs
buildSatDatabase()

# init the slack client
sc = SlackClient(slack_token)

# seo telescope
telescope = Telescope(False)

# main loop
while True:
    # connect to slack
    logme('Trying to connect to Slack...')
    if sc.rtm_connect():
        logme('Connected to Slack!')
        slack_connected = True
        # get list of users
        slack_users = get_users()
        # get a list of channels
        slack_channels = get_channels()
        # get the id of the slack channel assigned to telescope interface
        slack_channel = None
        for channel in slack_channels:
            if channel['name'] == slack_channel_name:
                slack_channel = channel['id']
        if slack_channel == None:
            abort('Error! Could not find #%s.' % slack_channel_name)
        # send welcome message
        doWelcome()
        logme('Listening for commands on #%s...' % slack_channel_name)
        # data loop
        while True:
            try:
                msgs = sc.rtm_read()  # returns array of json objects, e.g. return of json.loads()
            except:
                logme('Error! Connection with Slack was lost. Retrying...')
                sc.rtm_connect()
            # process incoming messages
            process_messages(msgs)
            # ping to ensure connection is intact
            if (datetime.datetime.now() - dt_last_activity).total_seconds() > 60:
                logme('Pinging Slack server...')
                if not ping():
                    logme('Error! Connection with Slack was lost. Retrying...')
                    sc.rtm_connect()
                logme('Received pong...still connected to Slack!')
                dt_last_activity = datetime.datetime.now()
            # wait
            time.sleep(read_delay_s)
    else:
        logme("Error! Slack connection failed. Retrying in %d seconds..." %
              (slack_token, reconnect_delay_s))
    time.sleep(reconnect_delay_s)
