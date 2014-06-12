
class Range(object):
    def __init__(self, min, max):
        self.min = Min(min,max)
        self.max = Max(min,max)

    def good(self, value):
        if self.min is not None and self.max is not None:
            if value < self.min: return False
            if value > self.max: return False
        elif self.min is None and self.max is not None:
            if value > self.max: return False
        elif self.min is not None and self.max is None:
            if value < self.min: return False
        elif self.min is None and self.max is None:
            pass
        return True

class AreaRange(object):
    def __init__(self, xmin, xmax, ymin, ymax):
        self.x = Range(xmin, ymax)
        self.y = Range(ymin, ymax)

    def good(self, x, y):
        return self.x.good(x) and self.y.good(y)

class RequestBuilder(object):
    def __init__(self):
        pass

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

    def stationBased(self, t0, t1,
                    targetSamplingRate,
                    allowedGainCodes,
                    networkStationCodes,
                    stationRestrictionArea,

                    eventRestrictionArea = None,
                    magnitudeRange = None,
                    depthRange = None,
                    azimutalRange = None
                    ):
        pass

    def eventBased(self, t0, t1,
                   targetSamplingRate,
                   allowedGainCodes,
                   eventRestrictionArea,
                   magnitudeRange,
                   depthRange,

                   networkStationCodes = None,
                   stationRestrictionArea = None,
                   azimutalRange = None
                   ):
        pass
