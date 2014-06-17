import sys
from obspy.fdsn import Client as fClient
from obspy.core import AttribDict
from obspy.core.utcdatetime import UTCDateTime
from obspy.core import util
from obspy import taup

class Range(object):
    def __init__(self, minvalue, maxvalue):
        self._min = min([minvalue,maxvalue])
        self._max = max([minvalue,maxvalue])

    def min(self):
        return self._min

    def max(self):
        return self._max

    def good(self, value):
        if self._min is not None and self._max is not None:
            if value < self._min: return False
            if value > self._max: return False
        elif self._min is None and self._max is not None:
            if value > self._max: return False
        elif self._min is not None and self._max is None:
            if value < self._min: return False
        elif self._min is None and self._max is None:
            pass
        return True

class AreaRange(object):
    def __init__(self, xmin, xmax, ymin, ymax):
        self.x = Range(xmin, xmax)
        self.y = Range(ymin, ymax)

    def xmin(self):
        return self.x.min()

    def xmax(self):
        return self.x.max()

    def ymin(self):
        return self.y.min()

    def ymax(self):
        return self.y.max()

    def good(self, x, y):
        return self.x.good(x) and self.y.good(y)

class RequestBuilder(object):
    def __init__(self, server):
        if isinstance(server, str):
            self.client = fClient(server)
        elif isinstance(server, fClient):
            self.client = server
        else:
            raise BadParameter("Invalid server object, expected String address of fdsnClient Class")

    '''
    Generic Methods that can be used in FDSN or ARCLINK modes
    '''

    def __find_arrivalTime(self, t0, delta, depth, phase = ["P", "Pn", "Pdiff", "PKP", "PKiKP", "PKPab", "PKPbc", "PKPdf", "PKPdiff"]):
        '''
            delta is in degrees
            depth is in meters
        '''
        times = taup.getTravelTimes(delta, depth / 1000.0)

        seltimes = filter(lambda x: x['phase_name'] in phase, times)
        if len(seltimes) == 0: 
            return (None, None)

        time = seltimes[0]['time']
        slowness = seltimes[0]['dT/dD']
        return (t0 + time, slowness)

    def __build_event_dictionary(self, t0, originLatitude, originLongitude, originDepth):
        '''
            t0 is a UTCDateTime
            eventLatitude is degrees
            eventLongitude is degrees
            eventDepth is in Meters
        '''
        eventinfo = {
                     'eventId': "%s" % t0.strftime("%y%m%d_%H%M%S"),
                     'time': t0,
                     'latitude': originLatitude,
                     'longitude': originLongitude,
                     'depth': originDepth / 1000.0
                     }
        return AttribDict(eventinfo)

    def __build_station_dictionary(self, networkCode, stationCode, stationLatitude, stationLongitude, stationElevation):
        '''
            networkCode is String
            stationCode is String
            stationLatitude is degrees
            stationLongitude is degrees
            stationElevation is Meters
        '''
        stationinfo = {
                       'stationId': "%s.%s" % (networkCode, stationCode),
                       'latitude' : stationLatitude,
                       'longitude': stationLongitude,
                       'elevation': stationElevation / 1000.0
                   }
        return AttribDict(stationinfo)

    def __build_pick_dictionary(self, phase, time, slowness):
        pickinfo = {
                    'phase': phase,
                    'time': time,
                    'slowness': slowness
                    }
        return AttribDict(pickinfo)

    def __organize_by_station(self, lines):
        request = {}

        for line in lines:
            key = "%s.%s" % (line[2],line[3])
            try:
                list = request[key]
            except KeyError:
                request[key] = []
                list = request[key]
            list.append(line)

        return request

    def __organize_by_event(self, lines):
        request= {}

        for line in lines:
            key = "%s" % (line[2])
            try:
                list = request[key]
            except KeyError:
                request[key] = []
                list = request[key]
            list.append(line)

        return request


    '''
    FDSN specific Methods
    '''

    def _chooseFDSN(self, item, channel, targetsps, instcode):

        # Build a newitem for the current channel
        newitem = (channel.location_code, channel.code, channel.sample_rate, channel.azimuth, channel.dip, abs(channel.sample_rate - targetsps))

        # Check that the instrument code is allowed to select
        if channel.code[1] not in instcode: return item

        # If we did not select yet return the current
        if item[0] == None: return newitem

        # Decided based on SPS first
        if newitem[5] == item[5]:
            # And based on the Channel Instrument Code
            if instcode.index(newitem[1][1]) < instcode.index(item[1][1]):
                return newitem
        elif newitem[5] < item[5] and newitem[5] > targetsps:
            return newitem

        return item

    def _getFDSNChannelList(self, station, t0, targetsps, instcode):
        clist = [ ]

        ## Choose will build something like: (loca, chan, sps , az  , dip , dsps)
        z = (None, None, None, None, None, None)
        n = (None, None, None, None, None, None)
        e = (None, None, None, None, None, None)

        for channel in station.channels:
            if t0 < channel.start_date or (channel.end_date is not None and t0 > channel.end_date):
                continue

            if channel.code[-1] == "Z":
                z = self._chooseFDSN(z, channel, targetsps, instcode)
            elif channel.code[-1] == "N" or channel.code[-1] == "1":
                n = self._chooseFDSN(n, channel, targetsps, instcode)
            elif channel.code[-1] == "E" or channel.code[-1] == "2":
                e = self._chooseFDSN(e, channel, targetsps, instcode)
#             else:
#                 print >>sys.stderr, "Unknow channel %s.%s" % (channel.location_code, channel.code)

        # Return only (location, channel, azimuth, dip)
        z = (z[0], z[1], z[3], z[4])
        n = (n[0], n[1], n[3], n[4])
        e = (e[0], e[1], e[3], e[4])

        return (z,n,e)


    '''
    Request Builder Methods
    '''

    def stationBased(self, t0, t1, targetSamplingRate, allowedGainCodes, timeRange,
                    networkStationCodes,
                    stationRestrictionArea,

                    eventRestrictionArea = None,
                    magnitudeRange = None,
                    depthRange = None,
                    distanceRange = None):

        # List of request lines
        lines = []

        # Start Build the parameters to pass on to the Client
        kwargsstation = {
                  "starttime": t0,
                  "endtime": t1,
                  "level": "channel"
                  }

        # Check if network/station code is a string -> list
        if isinstance(networkStationCodes, str):
            networkStationCodes = [networkStationCodes]

        # Add rectangular coordinate restrictions
        if stationRestrictionArea:
            kwargsstation["minlatitude"] = stationRestrictionArea.ymin()
            kwargsstation["maxlatitude"] = stationRestrictionArea.ymax()
            kwargsstation["minlongitude"] = stationRestrictionArea.xmin()
            kwargsstation["maxlongitude"] = stationRestrictionArea.xmax()

        ## Start the Loop
        # On the given station patterns
        for code in networkStationCodes:
            (net, sta) = code.split(".")
            kwargsstation['net'] = net
            kwargsstation['sta'] = sta
            inventory = self.client.get_stations(**kwargsstation)

            # On networks found
            for network in inventory.networks:
                # On Stations found
                for station in network.stations:
                    # Start Build the parameters to pass on the Event Client
                    kwargsevent = {
                               "starttime": max([station.start_date, t0]),
                               "endtime": min([station.end_date, t1])
                               }

                    # Add the event area restrictions
                    if eventRestrictionArea:
                        kwargsevent["minlatitude"] = eventRestrictionArea.ymin()
                        kwargsevent["maxlatitude"] = eventRestrictionArea.ymax()
                        kwargsevent["minlongitude"] = eventRestrictionArea.xmin()
                        kwargsevent["maxlongitude"] = eventRestrictionArea.xmax()
                    elif distanceRange:
                        kwargsevent["latitude"]  = station.latitude
                        kwargsevent["longitude"] = station.longitude
                        kwargsevent["minradius"] = distanceRange.min()
                        kwargsevent["maxradius"] = distanceRange.max()

                    # Event Depth Filter
                    if depthRange:
                        if depthRange.min() is not None:
                            kwargsevent["mindepth"] = depthRange.min()
                        if depthRange.max() is not None:
                            kwargsevent["maxdepth"] = depthRange.max()

                    # Magnitude Filter
                    if magnitudeRange:
                        if magnitudeRange.min() is not None:
                            kwargsevent["minmagnitude"] = magnitudeRange.min()
                        if magnitudeRange.max() is not None:
                            kwargsevent["maxmagnitude"] = magnitudeRange.max()

                    events = self.client.get_events(**kwargsevent)

                    i=0
                    # Event loop
                    if len(events) > 0:
                        for event in events:
                            i += 1
                            origin = event.preferred_origin()
                            if origin is None:
                                print >>sys.stderr,"Bad origin."
                                continue
                            delta = util.locations2degrees(station.latitude,
                                                           station.longitude,
                                                           origin.latitude,
                                                           origin.longitude)
                            (ta, slowness) = self.__find_arrivalTime(origin.time, delta, origin.depth)
                            (z,n,e) = self._getFDSNChannelList(station, ta, targetSamplingRate, allowedGainCodes)

                            EI = self.__build_event_dictionary(origin.time, origin.latitude, origin.longitude, origin.depth)
                            SI = self.__build_station_dictionary(network.code, station.code, station.latitude, station.longitude, station.elevation)
                            PI = self.__build_pick_dictionary("P", ta, slowness)
                            
                            lines.append((ta + timeRange.min(),
                                          ta + timeRange.max(),
                                          network.code,
                                          station.code,
                                          [z,n,e],
                                          SI,
                                          EI,
                                          PI)
                                         )
                    else:
                        print >>sys.stderr,"No Events Found"

        request = self.__organize_by_station(lines)

        for key in request:
            print key
            for line in lines:
                print " %s - %s\n %s %s\n %s\n %s\n %s\n %s" % line
                print ""
            print ""


    def eventBased(self, t0, t1, targetSamplingRate, allowedGainCodes, timeRange,
                   eventRestrictionArea,
                   magnitudeRange,
                   depthRange,

                   networkStationCodes = None,
                   stationRestrictionArea = None,
                   distanceRange = None
                   ):
        lines = []
        
        
        kwargsevent = {
                       "starttime": t0,
                       "endtime": t1,
                       }
        
        if eventRestrictionArea:
            kwargsevent["minlatitude"] = eventRestrictionArea.ymin()
            kwargsevent["maxlatitude"] = eventRestrictionArea.ymax()
            kwargsevent["minlongitude"] = eventRestrictionArea.xmin()
            kwargsevent["maxlongitude"] = eventRestrictionArea.xmax()

        if magnitudeRange:
            kwargsevent["minmagnitude"] = magnitudeRange.min()
            kwargsevent["maxmagnitude"] = magnitudeRange.max()

        if depthRange:
            kwargsevent["mindepth"] = depthRange.min()
            kwargsevent["maxdepth"] = depthRange.max()

        events = self.client.get_events(**kwargsevent)

        for event in events:
            origin = event.preferred_origin()

            if origin is None:
                print >>sys.stderr,"Bad origin."
                continue
            
            kwargsstation = {
                   "starttime": (origin.time - 86400),
                   "endtime": (origin.time +86400),
                   "level": "channel"
                   }
            
            if stationRestrictionArea:
                kwargsstation["minlatitude"] = stationRestrictionArea.ymin()
                kwargsstation["maxlatitude"] = stationRestrictionArea.ymax()
                kwargsstation["minlongitude"] = stationRestrictionArea.xmin()
                kwargsstation["maxlongitude"] = stationRestrictionArea.xmax()

            if distanceRange:
                kwargsstation["latitude"]  = origin.latitude
                kwargsstation["longitude"] = origin.longitude
                kwargsstation["minradius"] = distanceRange.min()
                kwargsstation["maxradius"] = distanceRange.max()


            if networkStationCodes:
                for code in networkStationCodes:
                    (net, sta) = code.split(".")
                    kwargsstation['net'] = net
                    kwargsstation['sta'] = sta
#                     print net,sta
                    inventory = self.client.get_stations(**kwargsstation)
                    i=0
                    # Event loop
                    if inventory is not None and len(inventory.networks) > 0:
                        for network in inventory.networks:
                            for station in network.stations:
                                i += 1
                                delta = util.locations2degrees(station.latitude,
                                                           station.longitude,
                                                           origin.latitude,
                                                           origin.longitude)
                                (ta, slowness) = self.__find_arrivalTime(origin.time, delta, origin.depth)
                                (z,n,e) = self._getFDSNChannelList(station, ta, targetSamplingRate, allowedGainCodes)
        
                                EI = self.__build_event_dictionary(origin.time, origin.latitude, origin.longitude, origin.depth)
                                SI = self.__build_station_dictionary(network.code, station.code, station.latitude, station.longitude, station.elevation)
                                PI = self.__build_pick_dictionary("P", ta, slowness)
        
                                lines.append((ta + timeRange.min(),
                                          ta + timeRange.max(),
                                          EI['eventId'],
                                          [z,n,e],
                                          SI,
                                          EI,
                                          PI)
                                          )


                    else:
                        print >>sys.stderr,"No Stations Found"
            else:
                print >>sys.stderr, "Erro, sem estacao especificada use *.* para pedir tudo."
                sys.exit()

        request = self.__organize_by_event(lines)
        for key in request:
            print key
            for line in lines:
                print " %s - %s\n %s %s\n %s\n %s\n %s\n" % line
                print ""
            print ""




if __name__ == "__main__":
    # Get an Instance of the class
    rb = RequestBuilder("IRIS")

    # Call the stationBased
    rb.eventBased(t0 = UTCDateTime("2011-01-01"),
                    t1 = UTCDateTime("2011-02-01"),
                    targetSamplingRate = 20.0,
                    allowedGainCodes = ["H", "L"],
                    timeRange = Range(-20.0, 50.0),

                    networkStationCodes = [ "TA.A*"],
                    stationRestrictionArea = AreaRange(-150.0, -90.0, 15.0, 35.0),

                    eventRestrictionArea = AreaRange(-75.0, -15.0, -35.0, -45.0),
                    magnitudeRange = Range(5.0, 7.0),
                    depthRange = Range(0.0, 400.0),
                    distanceRange = None
                    )
