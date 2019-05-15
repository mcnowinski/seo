import datetime
import subprocess
import sys
import urllib2
import json

###############################
#CHANGE THESE VALUES AS NEEDED#
###############################
# log file
log_fname = 'clouds.log'
# dat file
dat_fname = 'clouds.' + datetime.datetime.now().strftime("%Y%m%d_%H%M%S") + '.dat'
# wunderground token
wunderground_token = sys.argv[1]
wunderground_station = 'KCASONOM51'


def logme(msg):
    # open log file
    log = open(log_fname, 'a+')
    dt = datetime.datetime.now().strftime("%Y/%m/%d %H:%M:%S:\t")
    log.write((dt + msg + "\n").encode('utf8'))
    log.close()
    print (dt + msg).encode('utf8')
    return


def datme(msg):
    # open log file
    dat = open(dat_fname, 'a+')
    dt = datetime.datetime.now().strftime("%Y/%m/%d %H:%M:%S:\t")
    dat.write((msg + "\n").encode('utf8'))
    dat.close()
    #print (dt + msg).encode('utf8')
    return

# get weather from Wunderground


def getForecast(command, user):
    logme('Retrieving the hourly forecast from wunderground.com...')

    f = urllib2.urlopen('http://api.wunderground.com/api/%s/conditions/q/pws:%s.json' %
                        (wunderground_token, wunderground_station))
    json_string = f.read()
    parsed_json = json.loads(json_string)
    #print json_string
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


def getCloudCover(station):
    row = {}

    logme('Retrieving the current cloud cover for %s.' % station)

    #match = re.search('^\\\\(weather)\s(hourly)', command)
    # just in case wunderground is down...
    try:
        url = 'http://api.wunderground.com/api/%s/conditions/q/pws:%s.json' % (
            wunderground_token, station)
        f = urllib2.urlopen(url)
        json_string = f.read()
        parsed_json = json.loads(json_string)
        for field in conditions_current_observation_fields:
            if len(field) == 1:
                row[field[0]] = parsed_json['current_observation'][field[0]]
            elif len(field) == 2:
                row[field[0]+'_'+field[1]] = parsed_json['current_observation'][field[0]][field[1]]
        data.append(row)
    except:
        logme('Error. Could not get cloud cover for %s.' % station)
        return


def runSubprocess(command_array):
    # command array is array with command and all required parameters
    try:
        sp = subprocess.Popen(
            command_array, stderr=subprocess.PIPE, stdout=subprocess.PIPE)
        sp.wait()
        logme('Running subprocess ("%s" %s)...' %
              (' '.join(command_array), sp.pid))
        output, error = sp.communicate()
        return (output, error, sp.pid)
    except:
        logme('Error. Subprocess ("%s" %d) failed.' %
              (' '.join(command_array), sp.pid))
        return ('', '', 0)


stations = ['KCASONOM51','KCASONOM64', 'KCASONOM6', 'KCASONOM43']
conditions_current_observation_fields = [('station_id',), ('weather',),
                              ('precip_1hr_in',), ('precip_today_in',), ('temp_f',), ('dewpoint_f',), ('observation_location','latitude'), ('observation_location','longitude'), ('observation_location','elevation')]
data = []
for station in stations:
    getCloudCover(station)
for item in data:
    datme(','.join(map(str,item.values())))
