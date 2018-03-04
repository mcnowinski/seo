import os
import glob
import math
import subprocess
import re
import sys
import string
from decimal import Decimal
from astropy.io import fits
from astropy.wcs import WCS
from dateutil import parser
import matplotlib.pyplot as plt
import numpy as np
import numpy.ma as ma
from scipy.ndimage import median_filter
from pyds9 import DS9
import argparse
import pandas as pd

# logger


def logme(str):
    log.write(str + "\n")
    print str
    return


def exit():
    logme('Program execution halted.')
    log.close()
    os.sys.exit(1)

# run external process


def runSubprocess(command_array):
    # command array is array with command and all required parameters
    try:
        with open(os.devnull, 'w') as fp:
            sp = subprocess.Popen(command_array, stderr=fp, stdout=fp)
        # logme('Running subprocess ("%s" %s)...'%(' '.join(command_array), sp.pid))
        sp.wait()
        output, error = sp.communicate()
        return (output, error, sp.pid)
    except:
        logme('Error. Subprocess ("%s" %d) failed.' %
              (' '.join(command_array), sp.pid))
        return ('', '', 0)


#
# START SETTINGS
# MODIFY THESE FIELDS AS NEEDED!
#
# input path *with* ending forward slash
input_path = './'
# output path *with* ending forward slash
sex_output_path = './sex/'
# suffix for output files, if any...
sex_output_suffix = '.sex'
# log file name
log_fname = './log.sextractor.txt'
# path to sextractor executable and config files (incl. the filenames!)
sextractor_bin_fname = 'C:/owncloud/code/seo/sextractor/sextractor.exe'
sextractor_cfg_fname = 'C:/owncloud/code/seo/sextractor/seo.sex'
sextractor_param_fname = 'C:/owncloud/code/seo/sextractor/seo.param'
sextractor_filter_fname = 'C:/owncloud/code/seo/sextractor/seo.conv'
# tolerance for object matching
dRa = 0.00062
dDec = 0.00062
# target/comp list
stars_in_fname = './stars.in.txt'
stars_out_fname = './stars.out.csv'
# mask file that identifies bad pixels
bad_pixels_fname = './bad_pixels.txt'
cleaned_output_path = './cor/'
#
# END SETTINGS
#

# open log file
log = open(log_fname, 'a+')

# set up the command line argument parser
parser = argparse.ArgumentParser(
    description='Perform photometry using sextractor.')
parser.add_argument('--plot_apd', action='store_true',
                    help='Plot average object magnitude vs. aperture diameter for all images.')
parser.add_argument('--plot_ds9', action='store_true',
                    help='Plot apertures for each image using DS9.')
parser.add_argument("--apd", help="Perform analysis for this aperture diameter (pixels).",
                    type=int)
args = parser.parse_args()

# make sure input files and folder exist
inputs = [input_path, sextractor_bin_fname, sextractor_cfg_fname,
          sextractor_param_fname,  sextractor_filter_fname, stars_in_fname]
for input in inputs:
    if not os.path.exists(input_path):
        logme('Error. The file or path (%s) does not exist.' % input)
        exit()

# does output directory exist? If not, create it...
outputs = [sex_output_path, cleaned_output_path]
for output in outputs:
    try:
        os.mkdir(output)
    except:
        pass

# do we need to clean up any bad pixels?
bad_pixels = []
if os.path.exists(bad_pixels_fname):
    logme('Found list of bad pixels in %s.' % bad_pixels_fname)
    with open(bad_pixels_fname, 'r') as f:
        lines = f.readlines()
        for line in lines:
            xy = line.split()
            bad_pixels.append([int(xy[0]), int(xy[1])])
# if pixels have been marked as bad, replace their values with median values
if len(bad_pixels):
    # get a list of all FITS files in the input directory
    fits_files = glob.glob(input_path+'*.fits')+glob.glob(input_path+'*.fit')
    for fits_file in fits_files:
        # get data from fits file
        fits_data = fits.open(fits_file)
        data = fits_data[0].data
        # initialize mask with zeros
        mask = np.zeros((np.size(data, 0), np.size(data, 1)), dtype=bool)
        # mark bad pixels in mask
        for bad_pixel in bad_pixels:
            mask[bad_pixel[1]-1][bad_pixel[0]-1] = 1
        # calculate median filtered version of data (using surrounding pixels)
        data_cleaned = np.copy(data)
        data_cleaned[mask] = median_filter(
            data, footprint=[[1, 1, 1], [1, 0, 1], [1, 1, 1]])[mask]
        # replace fits data
        fits_data[0].data = data_cleaned
        # write it to cleaned fits file
        output_file = os.path.basename(fits_file)
        output_file = cleaned_output_path+output_file
        fits_data.writeto(output_file, overwrite=True)
    # use corrected images as input to sextractor
    input_path = cleaned_output_path

# grab aperture settings, hopefully it will match the magnitude array ;)
apertures = []
apertures_string = ''
with open(sextractor_cfg_fname) as f:
    lines = f.readlines()
    for line in lines:
        match = re.match(r'^PHOT_APERTURES([\s\.0-9\,]+)', line)
        if match:
            apertures_string = match.group(1).strip()
            apertures = np.array(
                [int(aperture) for aperture in apertures_string.split(',')])
logme('Photometry to be performed for %d aperture diameters: %s.' %
      (len(apertures), apertures_string))

image_data = []
# get a list of all FITS files in the input directory
fits_files = glob.glob(input_path+'*.fits')+glob.glob(input_path+'*.fit')
# loop through all qualifying files and perform sextraction
for fits_file in sorted(fits_files):
    fits_data = fits.open(fits_file)
    header = fits_data[0].header
    airmass = header['AIRMASS']
    try:
        JD = header['MJD-OBS']
    except KeyError:
        JD = header['JD']
    logme("Sextracting %s" % (fits_file))
    output_file = sex_output_path + \
        fits_file.replace('\\', '/').rsplit('/', 1)[1]
    output_file = '%s%s.txt' % (output_file, sex_output_suffix)
    # add input filename, output filename, airmass, and jd to sex_file list
    image_data.append(
        {'image': fits_file, 'sex': output_file, 'jd': JD, 'airmass': airmass})
    # sextract this file
    (output, error, id) = runSubprocess([sextractor_bin_fname, fits_file, '-c', sextractor_cfg_fname, '-catalog_name',
                                         output_file, '-parameters_name', sextractor_param_fname, '-filter_name', sextractor_filter_fname])
    if error:
        logme('Error. Sextractor failed: %s' % output)
        exit()
logme('Sextracted %d files.' % len(image_data))

# read ra/dec from target/comp stars list
object_data = []
sfile = file('%s' % stars_in_fname, 'rt')
lines = [s for s in sfile if len(s) > 2 and s[0] != '#']
sfile.close()
count = 0
target_index = -1
for index, l in enumerate(lines):
    spl = l.split()
    ra = float(spl[0])
    dec = float(spl[1])
    name = spl[2]
    if spl[3] == 'T' or spl[2] == 't':
        if target_index >= 0:
            logme('Error. More than one target was identified in %s!' %
                  stars_in_fname)
            os.sys.exit(1)
        target_index = index
    object_data.append(
        {'index': index, 'ra': ra, 'dec': dec, 'object_name': name})
if target_index < 0:
    logme('Error. No targets were specified in %s!' % stars_in_fname)
    os.sys.exit(1)
logme('Searching for %d objects in sextracted data.' % len(object_data))

# look for target/comp matches in sextracted files
sex_data = []
for image in image_data:
    found = 0
    lines = [s for s in file(image['sex'], 'rt') if len(s) > 2]
    for l in lines:
        spl = l.split()
        ra = float(spl[0])
        dec = float(spl[1])
        for s in object_data:
            if abs(ra-s['ra']) < dRa and abs(dec-s['dec']) < dDec:
                sex_data_element = {'object_index': s['index'], 'object_name': s['object_name'], 'object_ra': s['ra'], 'object_dec': s[
                    'dec'], 'jd': image['jd'], 'airmass': image['airmass'], 'image': image['image'], 'sex': image['sex']}
                sex_data_element['ra'] = spl[0]
                sex_data_element['dec'] = spl[1]
                sex_data_element['x'] = spl[2]
                sex_data_element['y'] = spl[3]
                sex_data_element['num_apertures'] = len(apertures)
                for i in range(0, len(apertures)):
                    sex_data_element['mag%02d' % apertures[i]] = spl[4+i]
                    sex_data_element['magerr%02d' %
                                     apertures[i]] = spl[4+len(apertures)+i]
                # print sex_data_element
                sex_data.append(sex_data_element)
                found += 1
    if found != len(object_data):
        logme('Warning! Found %s of %s objects in %s.' %
              (found, len(object_data), f[0]))
# print sex_data
logme('Found %d observations of %d objects in %d sextracted files.' %
      (len(sex_data), len(object_data), len(image_data)))

# save compiled sex data to a new file
ofile = file(stars_out_fname, 'wt')
line = 'index,name,airmass,jd'
jd_index = 3
mag_start_index = len(line.split(','))
for i in range(0, len(apertures)):
    line += ',mag%02d' % apertures[i]
for i in range(0, len(apertures)):
    line += ',magerr%02d' % apertures[i]
ofile.write(line+'\n')
# sort by star desig, then JD
sex_data = sorted(sex_data, key=lambda x: (x['object_name'], x['jd']))
for o in sex_data:
    line = '%d,%s,%f,%f' % (
        o['object_index'], o['object_name'], o['airmass'], o['jd'])
    for i in range(0, len(apertures)):
        line += ',%s' % o['mag%02d' % apertures[i]]
    for i in range(0, len(apertures)):
        line += ',%s' % o['magerr%02d' % apertures[i]]
    ofile.write(line+'\n')
ofile.close()

# plot average mag vs aperture diameter if requested
if args.plot_apd:
    ofile = file(stars_out_fname, 'r')
    data = np.genfromtxt(ofile, delimiter=',', skip_header=1)
    for index, s in enumerate(object_data):
        filtered_array = np.array(filter(lambda row: row[0] == index, data))
        magnitudes = np.mean(filtered_array, axis=0)[
            mag_start_index:mag_start_index+len(apertures)]
        magnitude_stdevs = np.std(filtered_array, axis=0)[
            mag_start_index:mag_start_index+len(apertures)]
        # apertures = np.linspace(1,len(magnitudes),len(magnitudes))
        # plt.errorbar(apertures, magnitudes, yerr=magnitude_stdevs, marker='o',
        #             color='black', elinewidth=0.5, linestyle='None', markersize=3)
        plt.errorbar(apertures, magnitudes, marker='o',
                     color='black', elinewidth=0.5, linestyle='None', markersize=3)
        plt.gca().invert_yaxis()
        # plt.gca().axes.xaxis.set_ticklabels([])
        plt.xlabel('Aperture Diameter, D (pixels)')
        plt.ylabel('Ave. Instrumental Magnitude, m')
        plt.title(s['object_name'])
        plt.show()

# plot target and comp stars in ds9
if args.plot_ds9:
    ds = DS9()
    ds.set('frame clear #all')
    for image in image_data:
        fits_file = image['image']
        fname = os.path.abspath(fits_file).replace('\\', '/')
        ds.set('file %s' % fname)
        ds.set('zscale')
        # ds.set('zoom to fit')
        w2 = WCS('%s' % fname)
        for s in object_data:
            # find a match in the sex data
            for o in sex_data:
                if o['image'] == fits_file and o['object_name'] == s['object_name']:
                    xp = int(float(o['x']))
                    yp = int(float(o['y']))
                    # print o['object_name'].replace('-','')
                    reg2 = 'regions command "point %s %s #color=lightgreen text=%s point=cross"' % (xp, yp, o['object_name'].replace(
                        '-', '').replace('.', ''))  # the times two because xvar is up and then again that value down
                    ds.set('%s' % (reg2))
        ds.set('frame new')
    ds.set('frame first')

if args.apd:
    logme('Analyzing photometry for aperture diameter = %d pixels.' % args.apd)
    # make sure this aperture is in our data set!
    apd_index = -1
    for index, aperture in enumerate(apertures):
        if aperture == int(args.apd):
            apd_index = mag_start_index + index
    if apd_index == -1:
        logme('Error! Could not find data for the aperture diameter = %s.' % args.apd)
        os.sys.exit(1)
    ofile = file(stars_out_fname, 'r')
    data = np.genfromtxt(ofile, delimiter=',', skip_header=1)
    for index, s in enumerate(object_data):
        filtered_array = np.array(filter(lambda row: row[0] == index, data))
        # print 'Target %s' % s['object_name']
        jds = filtered_array[:, jd_index]
        magnitudes = filtered_array[:, apd_index]
        magnitude_errors = filtered_array[:, apd_index+len(apertures)]
        plt.errorbar(jds, magnitudes, yerr=magnitude_errors, marker='o',
                     color='black', elinewidth=0.5, linestyle='None', markersize=3)
        plt.gca().invert_yaxis()
        plt.xlabel('Julian Date')
        plt.ylabel('Instrumental Magnitude, m')
        plt.title(s['object_name'])
        plt.show()

    # get average of comps
    filtered_array = np.array(filter(lambda row: row[0] != target_index, data))[
        :, [jd_index, mag_start_index+apd_index]]
    sorted_filtered_array = filtered_array[np.lexsort(
        np.transpose(filtered_array)[::-1])]
    #print sorted_filtered_array
    #print len(sorted_filtered_array)
    df = pd.DataFrame(sorted_filtered_array)
    ave_comp_magnitudes = df.groupby(
        np.arange(len(df))//(len(object_data)-1)).mean().values
    #print ave_comp_magnitudes
    #print len(ave_comp_magnitudes)
    # print np.mean(np.array(filtered_array).reshape(-1, len(object_data)-1), axis=1)
    plt.errorbar(ave_comp_magnitudes[:, 0], ave_comp_magnitudes[:, 1], marker='o',
                 color='black', elinewidth=0.5, linestyle='None', markersize=3)
    plt.gca().invert_yaxis()
    plt.xlabel('Julian Date')
    plt.ylabel('Instrumental Magnitude, m')
    plt.title('Average Comparison Star')
    plt.show()
