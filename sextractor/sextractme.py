import os
import glob
import math
import subprocess
import re
import sys
import string
from decimal import Decimal
from astropy.io import fits
from astropy import wcs
from dateutil import parser


def logme(str):
    log.write(str + "\n")
    print str
    return


def runSubprocess(command_array):
    # command array is array with command and all required parameters
    try:
        with open(os.devnull, 'w') as fp:
            sp = subprocess.Popen(command_array, stderr=fp, stdout=fp)
        #logme('Running subprocess ("%s" %s)...'%(' '.join(command_array), sp.pid))
        sp.wait()
        output, error = sp.communicate()
        return (output, error, sp.pid)
    except:
        logme('Error. Subprocess ("%s" %d) failed.' %
              (' '.join(command_array), sp.pid))
        return ('', '', 0)


# MODIFY THESE FIELDS AS NEEDED!
# input path *with* ending forward slash
input_path = './'
# output path *with* ending forward slash
sex_output_path = './sex/'
# suffix for output files, if any...
sex_output_suffix = '.sex'
# log file name
log_fname = './log.sextractor.txt'
# path to sextractor executable and config files (incl. the filenames!)
sextractor_bin_fname = 'C:/owncloud/code/variable/sextractor/sextractor.exe'
sextractor_cfg_fname = 'C:/owncloud/code/variable/sextractor/seo.sex'
sextractor_param_fname = 'C:/owncloud/code/variable/sextractor/seo.param'
sextractor_filter_fname = 'C:/owncloud/code/variable/sextractor/seo.conv'
# tolerance for object matching
dRa = 0.00062
dDec = 0.00062
# target/comp list
stars_in_fname = './stars.in.txt'
stars_out_fname = './stars.out.csv'

# make sure input files and folder exist
if not os.path.exists(input_path):
    logme('Error. The path %s does not exist.' % input_path)
    os.sys.exit(1)
if not os.path.exists(sextractor_bin_fname):
    logme('Error. The path %s does not exist.' % sextractor_bin_fname)
    os.sys.exit(1)
if not os.path.exists(sextractor_cfg_fname):
    logme('Error. The path %s does not exist.' % sextractor_cfg_fname)
    os.sys.exit(1)
if not os.path.exists(sextractor_param_fname):
    logme('Error. The path %s does not exist.' % sextractor_param_fname)
    os.sys.exit(1)
if not os.path.exists(sextractor_filter_fname):
    logme('Error. The path %s does not exist.' % sextractor_filter_fname)
    os.sys.exit(1)
if not os.path.exists(stars_in_fname):
    logme('Error. The path %s does not exist.' % stars_in_fname)
    os.sys.exit(1)

# do output directories exist? If not, create them...
try:
    os.mkdir(sex_output_path)
except:
    pass

# open log file
log = open(log_fname, 'a+')

sex_files = []
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
    #print airmass, JD
    logme("Sextracting %s" % (fits_file))
    output_file = sex_output_path + \
        fits_file.replace('\\', '/').rsplit('/', 1)[1]
    output_file = '%s%s.txt' % (output_file, sex_output_suffix)
    # add filename, airmass, and jd to sex_file list
    sex_files.append([output_file, airmass, JD])
    # sextract this file
    (output, error, id) = runSubprocess([sextractor_bin_fname, fits_file, '-c', sextractor_cfg_fname, '-catalog_name',
                                         output_file, '-parameters_name', sextractor_param_fname, '-filter_name', sextractor_filter_fname])
    if error:
        logme('Error. Sextractor failed: %s' % output)
        os.sys.exit(1)
logme('Sextracted %d files.' % len(sex_files))

# read ra/dec from target/comp stars list
sfile = file('%s' % stars_in_fname, 'rt')
lines = [s for s in sfile if len(s) > 2]
sfile.close()
starslist = []
for l in lines:
    spl = l.split()
    ra = float(spl[0])
    dec = float(spl[1])
    name = spl[2]
    starslist.append([ra, dec, name])
logme('Searching for %d objects...' % len(starslist))
#print starslist

# look for target/comp matches in sextracted files
outlist = []
for f in sex_files:
    found = 0
    lines = [s for s in file(f[0], 'rt') if len(s) > 2]
    for l in lines:
        spl = l.split()
        ra = float(spl[0])
        dec = float(spl[1])
        for s in starslist:
            if abs(ra-s[0]) < dRa and abs(dec-s[1]) < dDec:
                outlist.append('%s %s %s %s' % (s[2], f[1], f[2], l))
                found += 1
    if found != len(starslist):
        logme('Warning! Found %s of %s objects in %s.' %
              (found, len(starslist), f[0]))
#print outlist
logme('Found %d observations of %d objects in %d sextracted files.' %
      (len(outlist), len(starslist), len(sex_files)))

# save steller list in a new file
ofile = file(stars_out_fname, 'wt')
for s in outlist:
    s = re.sub(r'\s+', r',', s.strip())
    ofile.write(s+'\n')
ofile.close()
