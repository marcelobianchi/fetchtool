from obspy.fdsn import Client as fClient
from obspy.core import AttribDict
from obspy.core.utcdatetime import UTCDateTime

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

    def __find_arrivalTime(self, delta, depth, phase = ["P", "Pn", "Pdiff", "PKP", "PKiKP", "PKPab", "PKPbc", "PKPdf", "PKPdiff"]):
        '''
            delta is in degrees
            depth is in meters
        '''
#        origin = event.preferred_origin()
#        if origin is None:
#            origin = event.origins[0]
#        delta = util.locations2degrees(station.latitude, station.longitude, origin.latitude, origin.longitude)
        times = taup.getTravelTimes(delta, depth / 1000.0)

        seltimes = filter(lambda x: x['phase_name'] in phase, times)
        if len(seltimes) == 0: 
            return (None, None)

        time = seltimes[0]['time']
        return (origin.time + time, times[0]['dT/dD'])

    def __build_event_dictionary(self, t0, eventLatitude, eventLongitude, eventDepth):
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

    def stationBased(self, t0, t1, targetSamplingRate, allowedGainCodes,
                    networkStationCodes,
                    stationRestrictionArea,

                    eventRestrictionArea = None,
                    magnitudeRange = None,
                    depthRange = None,
                    distanceRange = None):

        # Start Build the parameters to pass on to the Client
        kwargsstation = {
                  "starttime": t0,
                  "endtime": t1,
                  "level": "station"
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
        for code in networkStationCodes:
            (net, sta) = code.split(".")
            kwargsstation['net'] = net
            kwargsstation['sta'] = sta
            inventory = self.client.get_stations(**kwargsstation)

            for network in inventory.networks:
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

                    print "Netowrk: ", network.code, " Station: ", station.code, " Found %d Events" % len(events)

    def eventBased(self, t0, t1, targetSamplingRate, allowedGainCodes,
                   eventRestrictionArea,
                   magnitudeRange,
                   depthRange,

                   networkStationCodes = None,
                   stationRestrictionArea = None,
                   distanceRange = None
                   ):
        pass



if __name__ == "__main__":
    # Get an Instance of the class
    rb = RequestBuilder("IRIS")

    # Call the stationBased
    rb.stationBased(t0 = UTCDateTime("2011-01-01"),
                    t1 = UTCDateTime("2011-02-01"),
                    targetSamplingRate = 20.0,
                    allowedGainCodes = ["H", "L"],

                    networkStationCodes = ["TA.A*", "TA.Z*"],
                    stationRestrictionArea = AreaRange(-150.0, -90.0, 15.0, 35.0),

                    eventRestrictionArea = AreaRange(-180.0, 0.0, 0.0, 90.0),
                    magnitudeRange = Range(6.0, 7.0),
                    depthRange = Range(0.0, 400.0),
                    distanceRange = None
                    )
