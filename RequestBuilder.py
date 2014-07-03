import sys
from obspy.fdsn import Client as fClient
from obspy.core import AttribDict
from obspy.core.utcdatetime import UTCDateTime
from obspy.core import util
from obspy import taup
from obspy.fdsn.header import FDSNException
import pickle
import os

class BadParameter(Exception):
    pass

class NextItem(Exception):
    pass

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

class BaseBuilder(object):

    '''
    Generic Methods that can be used in FDSN or ARCLINK modes
    '''

    def getChannelList(self, station, t0, targetsps, instcode):
        raise Exception("Implement your own override !")

    def build(self, lines, network, station, origin, phaselist, phasename, targetSamplingRate, allowedGainCodes, timeRange):
        delta = util.locations2degrees(station.latitude,
                                       station.longitude,
                                       origin.latitude,
                                       origin.longitude)

        (ta, slowness) = self.__find_arrivalTime(origin.time, delta, origin.depth, phaselist)

        # This method is overriden in each implementation
        (z,n,e) = self.getChannelList(station, ta, targetSamplingRate, allowedGainCodes)

        EI = self.__build_event_dictionary(origin.time, origin.latitude, origin.longitude, origin.depth)
        SI = self.__build_station_dictionary(network.code, station.code, station.latitude, station.longitude, station.elevation)
        PI = self.__build_pick_dictionary(phasename, ta, slowness)

        lines.append((ta + timeRange.min(),
                      ta + timeRange.max(),
                      network.code,
                      station.code,
                      [z,n,e],
                      SI,
                      EI,
                      PI)
                     )

        return

    def getOrigin(self, event):
        origin = event.preferred_origin()

        if origin == None:
            if len(event.origins) == 1:
                return event.origins[0]
            raise NextItem("Bad origin for event.")

        return origin

    def resolve_phasenames(self, value):
        plist = []

        '''
        Define more groups as needed:
        '''
        groups = {
          "pgroup": ["P", "Pn", "Pdiff", "PKP", "PKiKP", "PKPab", "PKPbc", "PKPdf", "PKPdiff" ]
        }

        if isinstance(value, str):
            value = [ value ]

        for item in value:
            try:
                plist.extend(groups[item])
            except KeyError:
                plist.append(item)

        return (plist[0], plist)

    def __find_arrivalTime(self, t0, delta, depth, phase):
        '''
            delta is in degrees
            depth is in meters
        '''
        times = taup.getTravelTimes(delta, depth / 1000.0)

        seltimes = filter(lambda x: x['phase_name'] in phase, times)
        if len(seltimes) == 0: 
            raise NextItem("No phases selected for Depth: %s Distance: %s Phases: %s" % (delta, depth, phase))

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

    def organize_by_station(self, lines):
        request = {}

        for line in lines:
            key = "%s.%s" % (line[2],line[3])
            try:
                clist = request[key]
            except KeyError:
                request[key] = []
                clist = request[key]
            clist.append(line)

        return request

    def organize_by_event(self, lines):
        request= {}

        for line in lines:
            key = "%s" % (line[6]['eventId'])
            try:
                clist = request[key]
            except KeyError:
                request[key] = []
                clist = request[key]
            clist.append(line)

        return request

    def fill_kwargsstation(self, t0, t1,
                             stationRestrictionArea,
                             latitude, longitude, distanceRange):

        kwargs = {
                  "starttime": t0,
                  "endtime": t1,
                  "level": "channel"
                  }

        # Add rectangular coordinate restrictions
        if stationRestrictionArea and distanceRange:
            raise Exception("Is not permitted to give parameters distanceRange and stationRestrictionArea")
        elif stationRestrictionArea is not None and distanceRange is None:
            kwargs["minlatitude"] = stationRestrictionArea.ymin()
            kwargs["maxlatitude"] = stationRestrictionArea.ymax()
            kwargs["minlongitude"] = stationRestrictionArea.xmin()
            kwargs["maxlongitude"] = stationRestrictionArea.xmax()
        elif distanceRange is not None and stationRestrictionArea is None:
            if latitude is None or longitude is None:
                raise Exception("Error setting the coordinate reference.")
            kwargs["latitude"]  = latitude
            kwargs["longitude"] = longitude
            kwargs["minradius"] = distanceRange.min()
            kwargs["maxradius"] = distanceRange.max()
        else:
            pass

        return kwargs

    def fill_kwargsevent(self, t0, t1,
                   eventRestrictionArea,
                   magnitudeRange,
                   depthRange,
                   latitude,
                   longitude,
                   distanceRange):

        kwargs = {
                       "starttime": t0,
                       "endtime": t1,
                       }
        if magnitudeRange:
            kwargs["minmagnitude"] = magnitudeRange.min()
            kwargs["maxmagnitude"] = magnitudeRange.max()

        if depthRange:
            kwargs["mindepth"] = depthRange.min()
            kwargs["maxdepth"] = depthRange.max()

        if distanceRange and eventRestrictionArea:
            raise Exception ("Is not permitted to give parameters distanceRange and eventRestrictionArea")
        elif distanceRange is not None and eventRestrictionArea is None:
            kwargs["latitude"]  = latitude
            kwargs["longitude"] = longitude
            kwargs["minradius"] = distanceRange.min()
            kwargs["maxradius"] = distanceRange.max()
        elif distanceRange is None and eventRestrictionArea is not None:
            kwargs["minlatitude"] = eventRestrictionArea.ymin()
            kwargs["maxlatitude"] = eventRestrictionArea.ymax()
            kwargs["minlongitude"] = eventRestrictionArea.xmin()
            kwargs["maxlongitude"] = eventRestrictionArea.xmax()
        else:
            pass

        return kwargs

    def show_request(self, request):
        for key in request:
            print key
            for line in request[key]:
                print " %s - %s\n %s %s\n %s\n %s\n %s\n %s\n" % line
                print ""
            print ""
        return

    def load_request(self, filename):
        if filename is None or not os.path.isfile(filename):
            raise Exception("Cannot read file, %s" % filename)

        try:
            iofile = file("request.pik", "r")
            request = pickle.load(iofile)
            iofile.close()
        except:
            request = None

        return request

    def save_request(self, filename, request):
        if filename is None or os.path.isfile(filename):
            raise Exception("Will not overwrite file, %s" % filename)

        if request is None: 
            raise Exception("Request cannot be empty")

        iofile = file(filename, "w")
        pickle.dump(request, iofile)
        iofile.close()

class RequestBuilder(BaseBuilder):
    def __init__(self, server):
        if isinstance(server, str):
            self.fdsn_client = fClient(server)
        elif isinstance(server, fClient):
            self.fdsn_client = server
        else:
            raise BadParameter("Invalid server object, expected String address of fdsnClient Class")

    '''
    FDSN specific Methods
    '''

    def __chooseFDSN(self, item, channel, targetsps, instcode):

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

    def getChannelList(self, station, t0, targetsps, instcode):

        ## Choose will build something like: (loca, chan, sps , az  , dip , dsps)
        z = (None, None, None, None, None, None)
        n = (None, None, None, None, None, None)
        e = (None, None, None, None, None, None)

        for channel in station.channels:
            if t0 < channel.start_date or (channel.end_date is not None and t0 > channel.end_date):
                continue

            if channel.code[-1] == "Z":
                z = self.__chooseFDSN(z, channel, targetsps, instcode)
            elif channel.code[-1] == "N" or channel.code[-1] == "1":
                n = self.__chooseFDSN(n, channel, targetsps, instcode)
            elif channel.code[-1] == "E" or channel.code[-1] == "2":
                e = self.__chooseFDSN(e, channel, targetsps, instcode)
#             else:
#                 print >>sys.stderr, "Unknow channel %s.%s" % (channel.location_code, channel.code)

        if z[0] is None or n[0] is None or e[0] is None:
            raise NextItem("No Z, N or E channel(s) found")

        # Return only (location, channel, azimuth, dip)
        z = (z[0], z[1], z[3], z[4])
        n = (n[0], n[1], n[3], n[4])
        e = (e[0], e[1], e[3], e[4])

        return (z,n,e)

    '''
    Request Builder Methods
    '''

    def stationBased(self, t0, t1, targetSamplingRate, allowedGainCodes, timeRange, phasesOrPhaseGroup, 
                    networkStationCodes,
                    stationRestrictionArea,

                    eventRestrictionArea = None,
                    magnitudeRange = None,
                    depthRange = None,
                    distanceRange = None):

        (phasename, phaselist) = self.resolve_phasenames(phasesOrPhaseGroup)
        print >>sys.stderr,"Searching using: %s %s" % (phasename, phaselist)

        # List of request lines
        lines = []

        # Start Build the parameters to pass on to the Client
        kwargsstation = self.fill_kwargsstation(t0,
                                                  t1,
                                                  stationRestrictionArea,
                                                  None, None, None)

        # Check if network/station code is a string -> list
        if isinstance(networkStationCodes, str):
            networkStationCodes = [ networkStationCodes ]

        ## Start the Loop
        # On the given station patterns
        for code in networkStationCodes:
            (net, sta) = code.split(".")
            kwargsstation['net'] = net
            kwargsstation['sta'] = sta

            try:
                inventory = self.fdsn_client.get_stations(**kwargsstation)
            except:
                inventory = None

            if inventory is None or len(inventory.networks) == 0:
                print >>sys.stderr,"No stations for pattern: %s.%s" % (net, sta)
                continue

            # On networks found
            print >>sys.stderr,"Stations for pattern: %s.%s" % (net,sta)
            for network in inventory.networks:
                # On Stations found
                for station in network.stations:
                    print >>sys.stderr,"\n Working on station %s.%s" % (network.code, station.code)

                    # Start Build the parameters to pass on the Event Client
                    kwargsevent = self.fill_kwargsevent(max([station.start_date, t0]),
                                                          min([station.end_date, t1]),
                                                          eventRestrictionArea,
                                                          magnitudeRange,
                                                          depthRange,
                                                          station.latitude,
                                                          station.longitude,
                                                          distanceRange)

                    try:
                        events = self.fdsn_client.get_events(**kwargsevent)
                    except FDSNException,e:
                        events = None

                    # Event loop
                    if events is None or len(events) < 0:
                        print >>sys.stderr,"  No Events Found."
                        continue

                    for event in events:
                        try:
                            origin = self.getOrigin(event)

                            print >>sys.stderr,"  Working on origin: %s" % str(origin.time),
                            self.build(lines,
                                       network, station, origin,
                                       phaselist, phasename,
                                       targetSamplingRate, allowedGainCodes, timeRange)
                            print >>sys.stderr,"OK!"
                        except NextItem,e:
                            print >>sys.stderr,"  Skipping: %s" % str(e)
                            continue

        request = self.organize_by_station(lines)

        return request

    def eventBased(self, t0, t1, targetSamplingRate, allowedGainCodes, timeRange, phasesOrPhaseGroup,
                   eventRestrictionArea,
                   magnitudeRange,
                   depthRange,

                   networkStationCodes = None,
                   stationRestrictionArea = None,
                   distanceRange = None
                   ):

        lines = []

        (phasename, phaselist) = self.resolve_phasenames(phasesOrPhaseGroup)
        print >>sys.stderr,"Searching using: %s %s" % (phasename, phaselist)

        if networkStationCodes is None:
            networkStationCodes = [ "*.*" ]

        kwargsevent = self.fill_kwargsevent(t0, t1,
                                              eventRestrictionArea,
                                              magnitudeRange,
                                              depthRange,
                                              None, None, None)

        try:
            events = self.fdsn_client.get_events(**kwargsevent)
            print >>sys.stderr,"Found %d events." % len(events)
        except FDSNException:
            print >>sys.stderr,"No events found for the given parameters."
            return None

        # Event loop
        for event in events:
            try:
                origin = self.getOrigin(event)
            except NextItem,e:
                print >>sys.stderr,"Skipping Origin: %s" % str(e)
                continue

            print >>sys.stderr,"Working on origin: %s" % str(origin.time)

            kwargsstation = self.fill_kwargsstation((origin.time - 86400),
                                                      (origin.time + 86400),
                                                      stationRestrictionArea,
                                                      origin.latitude,
                                                      origin.longitude,
                                                      distanceRange)

            for code in networkStationCodes:
                (net, sta) = code.split(".")
                kwargsstation['net'] = net
                kwargsstation['sta'] = sta

                try:
                    inventory = self.fdsn_client.get_stations(**kwargsstation)
                except FDSNException,e:
                    inventory = None

                if inventory is None or len(inventory.networks) == 0:
                    print >>sys.stderr,"\n No stations for pattern: %s.%s" % (net, sta)
                    continue

                print >>sys.stderr,"\n Stations for pattern: %s.%s" % (net,sta)
                # Event loop
                for network in inventory.networks:
                    for station in network.stations:
                        try:
                            print >>sys.stderr,"  Working on station %s.%s " % (network.code, station.code),
                            self.build(lines,
                                       network, station, origin,
                                       phaselist, phasename,
                                       targetSamplingRate, allowedGainCodes, timeRange)
                            print >>sys.stderr,"OK!"
                        except NextItem, e:
                            print >>sys.stderr,"\n  Skipping: %s" % str(e)

        request = self.organize_by_event(lines)

        return request


if __name__ == "__main__":
    # Get an Instance of the class
    rb = RequestBuilder("IRIS")

#     req = rb.load_request("request.pik")

    # Call the stationBased
    req = rb.eventBased(t0 = UTCDateTime("2007-01-01"),
                        t1 = UTCDateTime("2008-01-01"),
                        targetSamplingRate = 20.0,
                        allowedGainCodes = ["H", "L"],
                        timeRange = Range(-120, 600),
                        phasesOrPhaseGroup = "pgroup",

                        networkStationCodes = [ "TA.A*"],
                        stationRestrictionArea = AreaRange(-150.0, -90.0, 15.0, 60.0),

                        eventRestrictionArea = AreaRange(-35.0, -10.0, -55.0, -60.0),
                        magnitudeRange = Range(6.2, 9.0),
                        depthRange = Range(0.0, 400.0),
                        distanceRange = None
                        )

    rb.show_request(req)

#     rb.save_request("request.pik", req)
#     rb.show_request(req)

