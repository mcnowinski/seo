import os
import glob
import math
import subprocess
import re
import sys
from decimal import Decimal
import datetime
import getpass
from astropy.io import fits
from astropy import wcs
from astropy.io.fits import getheader
import astropy.coordinates as coord
import astropy.units as u

def logme( str ):
   log.write(datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S: ") + str + "\n")
   print str
   return
   
def quit( code ):
    #change filter back    
    logme('Changing filter *back* to %s.'%current_filter)
    os.system('pfilter %s'%current_filter)
    log.close()
    sys.exit(code)
   

##check command line parameters
#if len(sys.argv) > 3:
#    print 'usage: pinpoint <RA true (HH:MM:SS.SS)> <DEC true (DD:MM:SS.SS)>'
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
 
#print len(sys.argv)
 
#MODIFY THESE FIELDS AS NEEDED!
base_path='/tmp/'+datetime.datetime.now().strftime("%Y%m%d.%H%M%S%f.pinpoint.")
#log file name
log_fname = '/tmp/'+getpass.getuser()+'.pinpoint.log'
#path to astrometry.net solve_field executable
solve_field_path='/home/mcnowinski/astrometry/bin/solve-field'
#astrometry parameters
downsample=2
#bin=1, 0.75 arsec/pixel
scale_low=0.55
scale_high=2.00
radius=30.0 #up this to 30 deg, just in case scope is *way* off
cpu_limit=30
#offset limits (deg)
max_ra_offset=30.0
max_dec_offset=30.0
min_ra_offset=0.05
min_dec_offset=0.05
#how many pointing iterations to allow?
max_tries=5
#image command parameters
time=10
bin=2
fits_fname=base_path+'pointing.fits'

log=open(log_fname, 'a+')	 

current_filter = 'clear'
#ensure filter is clear!
#get current filter setting
output = subprocess.check_output('pfilter', shell=True) 
#match = re.search('([h\\-alpha|u\\-band||g\\-band|r\\-band|i\\-band|z\\-band|clear])', output)
match = re.search('([a-zA-Z0-1\\-]+)', output)
if match:
    current_filter = match.group(1)
    logme('Current filter is %s.'%current_filter)
    #set to clear (temporarily)
    logme('Changing filter setting to clear (temporarily).')    
    os.system('pfilter clear')
else:
    logme('Error. Unrecognized filter (%s).'%output)
 
ra_offset = 5.0
dec_offset = 5.0
iteration = 0
while((abs(ra_offset) > min_ra_offset or abs(dec_offset) > min_dec_offset) and iteration < max_tries):
    iteration += 1
    
    logme('Performing adjustment #%d...'%(iteration))

    #get pointing image  
    os.system('image time=%d bin=%d outfile=%s %s' % (time, bin, fits_fname, image_parms));
        
    if not os.path.isfile(fits_fname):
        logme('Error. File (%s) not found.' %fits_fname)
        quit(1)
 
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
            quit(1)
            
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
        quit(1)

    ra_offset=float(ra_target)-float(RA_image)
    if ra_offset > 350:
        ra_offset -= 360.0
    dec_offset=float(dec_target)-float(DEC_image)

    if(abs(ra_offset) <= max_ra_offset and abs(dec_offset) <=max_dec_offset):
        #os.system('tx offset ra=%f dec=%f cos > /dev/null' % (ra_offset, dec_offset))
        os.system('tx offset ra=%f dec=%f > /dev/null' % (ra_offset, dec_offset))
        logme("...complete (dRA=%f deg, dDEC=%f deg)."%(ra_offset, dec_offset))
    else:
        logme("Error. Calculated offsets too large (tx offset ra=%f dec=%f)! Pinpoint aborted." % (ra_offset, dec_offset))   
        quit(1)

if(iteration < max_tries):   
    logme('BAM! Your target has been pinpoint-ed!')
    quit(0)
    
logme('Error. Exceeded maximum number of adjustments (%d).'%max_tries)
quit(1)   