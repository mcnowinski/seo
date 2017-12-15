import os
import glob
import datetime
import pytz
import re
from collections import defaultdict
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.dates as dates
from matplotlib.dates import date2num
from dateutil import parser
from terminaltables import AsciiTable
from astropy.units import ccd
from statsmodels.discrete.tests.test_sandwich_cov import exposure
from matplotlib._path import points_in_path
from astropy import units as u
from astropy.coordinates import SkyCoord

#track weather
global weather
weather = []
#track users
global users
users = []
#track imaging
global images
images = []
#track pointing
global points
points = []
#start and end dates
global dtStart
global dtEnd

#
#this script parses aster log files, e.g. tserver.log
#

#regular expressions to identify aster log messages
re_doCCD = '^ccd time=(\S+) nbytes=[0-9]+ bin=([0-9]+) bzero=[0-9\\-]+\s?(dark)?'
re_doLock = '^done lock user=(\S+) email=(\S+)'
re_doSlit = '^slit (open|close)'
re_doPoint = '^point ra=(\S+) dec=(\S+)'
re_doWeather = '^done taux ovolts=[\S]+ irvolts=[\S]+ cloud=([0-9\\-\\+\\.]+) rain=(0|1) dew=(\S+)'
re_doWhere = '^done where ra=([\S]+) dec=([\S]+) equinox=[\S]+ ha=([\S]+) secz=([\S]+) alt=([\S]+) az=([\S]+)'

def calcPointing():
    #print points
    skycoords = SkyCoord(points, unit=(u.hourangle, u.deg))
    #print c
    ra_rad = []
    dec_rad = []
    for skycoord in skycoords:
        ra_rad.append(skycoord.ra.wrap_at(180 * u.deg).radian)
        dec_rad.append(skycoord.dec.radian)
    plt.figure(figsize=(8,4.2))
    plt.subplot(111, projection="aitoff")
    plt.title("SEO Targets from %s to %s"%(dtStart.strftime('%m-%d-%Y'), dtEnd.strftime('%m-%d-%Y')), y = 1.08)
    plt.grid(True)
    plt.plot(ra_rad, dec_rad, 'or', markersize=3)
    plt.subplots_adjust(top=0.95,bottom=0.0)
    #plt.show()
    plt.savefig('pointing.png', bbox_inches='tight')
    
def calcCCD():
    exposure_sum = exposure_light_sum = exposure_dark_sum = 0
    exposure_light_count = exposure_dark_count = 0
    for image in images:
        exposure = image[1]
        exposure_sum += exposure
        #binning = image[2]
        dark = image[3]
        if dark:
            exposure_dark_count += 1
            exposure_dark_sum += exposure
        else:
            exposure_light_count += 1
            exposure_light_sum += exposure
    logme('Total number of images taken was %d.'%len(images))
    logme(' * Number of light images taken was %d (%d%%).'%(exposure_light_count, int(float(exposure_light_count)/len(images)*100)))
    logme(' * Number of dark images taken was %d (%d%%).'%(exposure_dark_count, int(float(exposure_dark_count)/len(images)*100)))
    logme('Total imaging time was %d s.'%exposure_sum)  
    logme(' * Light imaging time was %d s (%d%%).'%(exposure_light_sum, int(float(exposure_light_sum)/exposure_sum*100))) 
    logme(' * Dark imaging time was %d s (%d%%).'%(exposure_dark_sum, int(float(exposure_dark_sum)/exposure_sum*100)))     

def calcUsage():
    tabledata = [] 
    tabledata.append(['User', 'Days Used', 'Last Used'])   
    
    for user in users:
        #for each user, get number of days online and last date online   
        uniqueDays = []
        dtLast = None
        for dt in user[2]:
            #awkward way to count all unique days this user locked the scope
            if len(uniqueDays) == 0:
                uniqueDays.append(dt.strftime("%m-%d-%Y"))
            else:
                foundDay = False
                for uniqueDay in uniqueDays:
                    if dt.strftime("%m-%d-%Y") == uniqueDay:
                        foundDay = True #already counted this day
                        break
                #add this day to the list    
                if not foundDay:
                    uniqueDays.append(dt.strftime("%m-%d-%Y"))
            #more awkwardness to track last date that this user was on            
            if dtLast == None:
                dtLast = dt
            elif dt > dtLast:
                dtLast = dt
        #print user[0], dtLast.strftime("%m-%d-%Y"), len(uniqueDays)
        tabledata.append([user[0], '%d'%len(uniqueDays), dtLast.strftime("%m-%d-%Y")])        
    table = AsciiTable(tabledata) 
    logme(table.table)

def calcWeather():
    #calc average cloud cover during observations
    cloud_average = 0
    uniqueDays = []
    for wx in weather:
        dt = wx[0]
        cloud = wx[1]
        rain = int(wx[2])
        #calc average cloud cover
        cloud_average += float(cloud)/len(weather)
        #calculate number of days with rain
        if rain == 1:
            if len(uniqueDays) == 0:
                uniqueDays.append(dt.strftime("%m-%d-%Y"))
            else:
                foundDay = False
                for uniqueDay in uniqueDays:
                    if dt.strftime("%m-%d-%Y") == uniqueDay:
                        foundDay = True #already counted this day
                        break
                #add this day to the list    
                if not foundDay:
                    uniqueDays.append(dt.strftime("%m-%d-%Y"))

    logme('Average cloud cover was %d%%.'%(int(cloud_average)))
    logme('Number of days with rain was %d (%d%%).'%(len(uniqueDays), int(float(len(uniqueDays))/(dtEnd-dtStart).days*100)))
    
    #plot cloud cover
    x = [date2num(date) for (date, cloud, rain) in weather]
    y = [cloud for (date, cloud, rain) in weather]
    fig = plt.figure()
    graph = fig.add_subplot(111)
    # Plot the data as a red line with round markers
    #graph.plot(x,y,'bo',markersize=1)
    graph.bar(x, y, width=0.05, color='b')
    # Set the xtick locations to correspond to just the dates you entered.
    graph.set_xticks(x)
    # Set the xtick labels to correspond to just the dates you entered.
    graph.set_xticklabels(
            [date.strftime("%m-%d") for (date, cloud, rain) in weather]
            )
    x_major_lct = dates.AutoDateLocator(minticks=2, maxticks=10, interval_multiples=True)
    #x_minor_lct = dates.HourLocator(byhour = range(0,25,1))
    x_fmt = dates.AutoDateFormatter(x_major_lct, defaultfmt='%m-%d')
    graph.xaxis.set_major_locator(x_major_lct)
    #graph.xaxis.set_minor_locator(x_minor_lct)
    graph.xaxis.set_major_formatter(x_fmt)
    plt.xlabel('Day')
    plt.ylabel('Cloud Cover (%)')
    #plt.show()
    plt.savefig('cloud.png', bbox_inches='tight')

def doWhere(msg, dt):
    match = re.search(re_doWhere, msg)
    if match:
        ra=match.group(1)
        dec=match.group(2)
        ha=match.group(3)
        secz=match.group(4)
        alt=match.group(5)
        az=match.group(6)
        #print 'Telescope was pointed to RA=%s, DEC=%s, ha=%s, secz=%s, alt=%s, az=%s.'%(match.group(1), match.group(2), match.group(3), match.group(4), match.group(5), match.group(6)) 
 
#re_doWeather = '^done taux ovolts=[\S]+ irvolts=[\S]+ cloud=([0-9\\-\\+\\.]+) rain=(0|1) dew=(\S+)'        
def doWeather(msg, dt):
    match = re.search(re_doWeather, msg)
    if match:
        #print match.group(1)
        cloud = int(float(match.group(1))*100)
        rain = match.group(2)
        #set negative cloudiness to zero
        if cloud < 0:
            cloud = 0
        #print 'Clouds=%s, rain=%s, dew=%s.'%(match.group(1), match.group(2), match.group(3))
        weather.append((dt, cloud, rain))         

def doPoint(msg, dt):
    match = re.search(re_doPoint, msg)
    if match:
        ra = match.group(1)
        dec = match.group(2)
        raDec = '%s %s'%(ra, dec)
        #make sure ra and dec values are valid; better way? regex?
        try:
            sc = SkyCoord(raDec, unit=(u.hourangle, u.deg))
            points.append(raDec)
        except:
            #logme('Error. Invalid RA or DEC value (%s).'%raDec)
            pass

def doSlit(msg, dt):
    match = re.search(re_doSlit, msg)
    if match:
        slit_status = match.group(1)
        #if match.group(1) == 'open':
        #    print 'Slit was opened at %s.'%dt.strftime("%Y-%m-%d %H:%M:%S")
        #else:
        #    print 'Slit was closed at %s.'%dt.strftime("%Y-%m-%d %H:%M:%S")            

def doLock(msg, dt):
    match = re.search(re_doLock, msg)
    if match:
        user=match.group(1)
        email=match.group(2)
        foundUser = False
        for user in users:
            if user[0] == email:
                foundUser = True
                user[2].append(dt)
                continue
        if not foundUser:
            users.append((email, user, [dt]))
        #print 'User=%s. Email=%s.'%(match.group(1), match.group(2))                

#re_doCCD = '^ccd time=(\S+) nbytes=[0-9]+ bin=([0-9]+) bzero=[0-9\\-]+\s?(dark)?'
def doCCD(msg, dt):
    match = re.search(re_doCCD, msg)
    if match:
        exposure=float(match.group(1))
        binning=int(match.group(2))
        if match.group(3):
            dark = 1
        else:
            dark = 0
        images.append((dt, exposure, binning, dark))
    #if match:
    #    #print msg
    #    if match.group(3):
    #        print 'Exposure time =\t%s sec. Binning = %s. Dark frame.'%(match.group(1), match.group(2))              
    #    else:
    #        print 'Exposure time =\t%s sec. Binning = %s. Light frame.'%(match.group(1), match.group(2))

def parse_message(msg, dt):
    match = False
    for msg_type in msg_types:
        #match the command first, then look for
        match = re.search(msg_type[0], msg)
        if match: #compare with list of known msgs
            #call associated function
            msg_type[1](msg, dt)
            break;
    #if not match: #did not recognize this command
    #    logme('Warning. Unrecognized message (%s).'%msg)

#valid message types
#element 1 is a regular expression describing the message (and its parameters)
#element 2 is the function called when this message is found
msg_types = [

#ccd time=2.310 nbytes=2097152 bin=2 bzero=32768
[re_doCCD, doCCD],
#done lock user=mcnowinski email=mcnowinski@gmail.com phone=7032869140 comment=fatflats timestamp=2017-02-28T02:02:23Z
[re_doLock, doLock],
#slit close
#slit open
[re_doSlit, doSlit],
#point ra=3:19:48.16 dec=41:30:42.103 equinox=2000.0
[re_doPoint, doPoint],
#done taux ovolts=2.867 irvolts=0.150 cloud=0.25 rain=0 dew=2.89
[re_doWeather, doWeather],
#done where ra=16:48:02.19 dec=-01:53:50.7 equinox=2017.166 ha=-35.853 secz=1.63 alt=38.0 az=132.0 slewing=0
[re_doWhere, doWhere],

]

#log the unlogger
def logme( str ):
   log.write(str + "\n")
   print str
   return 

#path to log files
input_path = '/data/'
#mask for log files
input_mask = 'tserver.*'
#log file name
log_fname = 'log.unlogger.txt'

log=open(log_fname, 'w')

#grab start and end dates if they are provided
#otherwise use look at 7 days
dtEnd = datetime.datetime.utcnow().replace(hour=0,minute=0,second=0)
dtStart = (dtEnd-datetime.timedelta(days=7)).replace(hour=23,minute=59,second=59)
#dtFilterStart = dtFilterEnd = None
#filterDates = False
if(len(os.sys.argv) >= 3):
    dtStart = parser.parse(os.sys.argv[1]).replace(hour=0,minute=0,second=0)
    dtEnd = parser.parse(os.sys.argv[2]).replace(hour=23,minute=59,second=59)
    #filterDates = True   
   
log_files=glob.glob(input_path+input_mask)

#dtStart = dtFilterStart
#dtEnd = dtFilterEnd
for log_file in log_files:
    #skip files that are obviously earlier than the date range we are interested in
    #tserver.2017-03-06T19:59:49Z
    match = re.search('([0-9]{4}\\-[0-9]{2}\\-[0-9]{2}T[0-9]{2}\\:[0-9]{2}\\:[0-9]{2}Z)', log_file) 
    if match:
        dt = parser.parse(match.group(1))
        if dt < pytz.utc.localize(dtStart): #this log file is before our time
            continue #go to next file
    #open file and read contents 
    #logme('Processing %s...'%log_file)    
    try:
        lf=open(log_file, 'r')
    except:
        #logme('Error. Failed to open log file (%s).'%log_file)
        continue      
    lines=lf.readlines()
    
    for line in lines:
        #logme(line)
        #grab only lines with a timestamp (for now)
        match = re.search('^([0-9]{4}\\-[0-9]{2}\\-[0-9]{2}T[0-9]{2}\\:[0-9]{2}\\:[0-9]{2}Z)', line) 
        #if it's a match, grab the datetimestamp and the log message   
        if match:
            #logme(match.group(1))
            dt = parser.parse(match.group(1))
            #print dt.strftime("%Y-%m-%d %H:%M:%S")
            #skip dates that have been filtered
            #if filterDates == True:
            if dt > pytz.utc.localize(dtEnd) or dt < pytz.utc.localize(dtStart):
                continue
            #else: #track start and end dates
            #    if dtStart == None:
            #        dtStart = dtEnd = dt #init start and end
            #    else:
            #        if dt < dtStart:
            #            dtStart = dt
            #        elif dt > dtEnd:
            #            dtEnd = dt
            msg = line[len(match.group(1)):len(line)].strip()
            #logme('%s\t%s'%(dt.strftime("%Y-%m-%d %H:%M:%S"),msg))
            #now let's parse the msg!
            #
            parse_message(msg, dt)
                
    lf.close()

logme('------------------------------------')
logme('Stone Edge Observatory (SEO) Metrics')
logme('------------------------------------')
logme('Date range is %s to %s.'%(dtStart.strftime("%m-%d-%Y"), dtEnd.strftime("%m-%d-%Y")))
  
#how was SEO weather?
calcWeather()

#what observations were performed
calcCCD()

#who used the scope and when
calcUsage() 

#where was the scope pointed?
calcPointing()

#close log
log.close()



