import os
import glob
import math
import subprocess
import re
import sys
from decimal import Decimal
import datetime
import getpass
##gain access to local astropy module
#sys.path.append('/home/mcnowinski/astropy/lib/python2.7/site-packages/astropy-1.2.1-py2.7-linux-x86_64.egg')
from astropy.io import fits
from astropy import wcs
from astropy.io.fits import getheader
import astropy.coordinates as coord
import astropy.units as u

def logme( str ):
   log.write(datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S: ") + str + "\n")
   print str
   return

##check command line parameters
#if len(sys.argv) > 3:
#    print 'usage: fixpoint <RA true (HH:MM:SS.SS)> <DEC true (DD:MM:SS.SS)>'
#    sys.exit(1)

#get target ra and dec from command line, if available
ra_target = None
dec_target = None
if len(sys.argv) >= 3:
    #convert to decimal degrees
    ra_target = coord.Angle(sys.argv[1], unit=u.hour).degree
    dec_target = coord.Angle(sys.argv[2], unit=u.deg).degree
    
image_parms = ''
if len(sys.argv) >= 4:
    #get any parameters passed for image command, e.g. notel
    image_parms = sys.argv[3]
    
#MODIFY THESE FIELDS AS NEEDED!
base_path='/tmp/'+datetime.datetime.now().strftime("%Y%m%d.%H%M%S%f.fixpoint.")
#log file name
log_fname = '/tmp/'+getpass.getuser()+'.fixpoint.log'
#path to astrometry.net solve_field executable
solve_field_path='/home/mcnowinski/astrometry/bin/solve-field'
#astrometry parameters
downsample=2
scale_low=0.55
scale_high=2.00
radius=10.0 #up this to 5 deg, just in case scope is way off
cpu_limit=30
#offset limits (deg)
max_ra_offset=10.0
max_dec_offset=10.0
min_ra_offset=0.05
min_dec_offset=0.05
#how many pointing iterations to allow?
max_tries=20
#image command parameters
time=10
bin=2
fits_fname=base_path+'pointing.fits'

log=open(log_fname, 'a+')	

#fits_fname = sys.argv[1]   

ra_offset = 5.0
dec_offset = 5.0
    
#logme('Performing adjustment #%d...'%(iteration))

#reset
os.system('tx offset ra=%f dec=%f > /dev/null' % (0.0, 0.0))
os.system('tx zero last > /dev/null')

#get pointing image
#ensure filter is clear!
#get current filter setting
output = subprocess.check_output('pfilter', shell=True) 
#match = re.search('([h\\-alpha|u\\-band||g\\-band|r\\-band|i\\-band|z\\-band|clear])', output)
match = re.search('([a-zA-Z0-1\\-]+)', output)
if match:
    logme('Current filter is %s.'%match.group(1))
    #set to clear (temporarily)
    logme('Changing filter setting to clear (temporarily).')    
    os.system('pfilter clear')
else:
    logme('Error. Unrecognized filter (%s).'%output)

os.system('image time=%d bin=%d outfile=%s %s' % (time, bin, fits_fname, image_parms));

#change filter back    
if match:
    logme('Changing filter *back* to %s.'%match.group(1))
    os.system('pfilter %s'%match.group(1))
    
if not os.path.isfile(fits_fname):
    logme('Error. File (%s) not found.' %fits_fname)
    log.close()
    sys.exit(1)

#get FITS header, pull RA and DEC for cueing the plate solving
if(ra_target == None or dec_target == None):
    header = getheader(fits_fname)
    try:
        ra_target = header['RA']
        ra_target = coord.Angle(ra_target, unit=u.hour).degree # create an Angle object
        dec_target = header['DEC']
        dec_target = coord.Angle(dec_target, unit=u.deg).degree
    except KeyError:
        logme("Error. RA/DEC not found in input FITS header (%s)." % fits_fname)
        log.close()
        os.remove(fits_fname)
        sys.exit(1)
        
#plate solve this image, using RA/DEC from FITS header
output = subprocess.check_output(solve_field_path + ' --no-verify --overwrite --no-remove-lines --downsample %d --scale-units arcsecperpix --scale-low %f --scale-high %f --ra %s --dec %s --radius %f --cpulimit %d --no-plots '%(downsample,scale_low,scale_high,ra_target,dec_target,radius,cpu_limit)+fits_fname, shell=True)
#output = subprocess.check_output(solve_field_path + ' --no-fits2fits --overwrite --downsample 2 --guess-scale --ra %s --dec %s --radius 1.0 --cpulimit 30 --no-plots '%(ra,dec)+'%s'%(new), shell=True)
#output = subprocess.check_output(solve_field_path + ' --no-verify --overwrite --no-remove-lines --downsample %d --guess-scale --ra %s --dec %s --radius %f --cpulimit %d --no-plots '%(downsample,ra_target,dec_target,radius,cpu_limit)+fits_fname, shell=True)
#let's go simpler after issues with gross telescope mis-pointing 021417
#output = subprocess.check_output(solve_field_path + ' --no-verify --overwrite --no-remove-lines --downsample %d --cpulimit %d --no-plots '%(downsample,cpu_limit)+fits_fname, shell=True)

#output = subprocess.check_output(solve_field_path + ' --no-verify --overwrite --downsample %d --cpulimit %d --no-plots '%(downsample,cpu_limit)+fits_fname, shell=True)
log.write(output)

#remove astrometry.net temporary files
os.remove(fits_fname)
os.remove(base_path+'pointing-indx.xyls')
os.remove(base_path+'pointing.axy')
os.remove(base_path+'pointing.corr')	
os.remove(base_path+'pointing.match')
os.remove(base_path+'pointing.rdls')
os.remove(base_path+'pointing.solved')
os.remove(base_path+'pointing.wcs')
os.remove(base_path+'pointing.new')

#look for field center in solve-field output
match = re.search('Field center\: \(RA,Dec\) \= \(([0-9\-\.\s]+)\,([0-9\-\.\s]+)\) deg\.', output)
if match:
    RA_image=match.group(1).strip()
    DEC_image=match.group(2).strip()		
else:
    logme("Error. Field center RA/DEC not found in solve-field output!")
    sys.exit(1)

ra_offset=float(ra_target)-float(RA_image)
if ra_offset > 350:
    ra_offset -= 360.0
dec_offset=float(dec_target)-float(DEC_image)

if(abs(ra_offset) <= max_ra_offset and abs(dec_offset) <=max_dec_offset):
    #print 'Calculated offset is dRA=%f dDEC=%f. Apply this as default (y/n)?'%(ra_offset, dec_offset)
    #yesno = raw_input('')
    #if yesno == 'y' or yesno == 'Y':
    os.system('tx offset ra=%f dec=%f > /dev/null' % (ra_offset, dec_offset))
    os.system('tx zero last > /dev/null')
    logme("Default offset set (dRA=%f deg, dDEC=%f deg)!"%(ra_offset, dec_offset))
    #else:
    #    logme("No offset applied.")        
else:
    logme("Error. Calculated offsets too large (tx offset ra=%f dec=%f)! Pinpoint aborted." % (ra_offset, dec_offset))   
    sys.exit(1)
    
sys.exit(0)

    