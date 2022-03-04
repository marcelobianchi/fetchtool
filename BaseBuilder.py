'''
FetchTool Interactive, Mutli-module seismological mass downloader package.
Copyright (C) 2015  Marcelo Bianchi <m.bianchi@iag.usp.br>

This file is part of FetchTool.

FetchTool is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
'''

from __future__ import division, print_function

import sys

from obspy.core import AttribDict
from obspy.core.utcdatetime import UTCDateTime
from obspy import geodetics
from obspy import taup
from obspy.core import event as ObsEvent
import numpy as np

import pickle
import os

STATUS = AttribDict({
                     "unset": "unset",

                     "requested": "requested",
                     "downloaded": "downloaded",
                     "saved": "saved",

                     "temporary_error": "temporary_error",
                     "permanent_error": "permanent_error",
                     })

MINUTE = 60
HOUR   = 60 * MINUTE
DAY    = 24 * HOUR

class Status(object):
    def __init__(self):
        self.level = STATUS.unset
        self.comment = None

    def __iter__(self):
        return self

    def __len__(self):
        return 0

    def next(self):
        raise StopIteration

    def show(self):
        print("Status Level: ",self.level, "[",self.comment,"]")


class BadParameter(Exception):
    pass


class NextItem(Exception):
    pass


class Range(object):
    def __init__(self, minvalue, maxvalue):
        self._min = minvalue # min([minvalue,maxvalue]) if minvalue != None and maxvalue != None else minvalue
        self._max = maxvalue # max([minvalue,maxvalue]) if minvalue != None and maxvalue != None else maxvalue

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

    @staticmethod
    def ALLDISTS():
        return Range(0.0, 180.0)
    
    @staticmethod
    def ALLMAGS():
        return Range(0.0, 10.0)

    @staticmethod
    def ALLDEPTHS():
        return Range(0.0, 1000.0)

    def __str__(self):
        return f"Range from: {self.min()} to: {self.max()}"


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

    @staticmethod
    def WORLD():
        return AreaRange(-180., 180., -90., 90)

    def __str__(self):
        return f'AreaRange from {self.xmin()} to {self.xmax()} and {self.ymin()} to {self.ymax()}'


class BaseBuilder(object):
    def __init__(self):
        self._plotevents = False
        self.__tworker = taup.TauPyModel()

    '''
    Generic Methods that can be used in FDSN or ARCLINK modes
    '''

    @staticmethod
    def _getChannelList(station, t0, targetsps, instcode):
        raise Exception("Implement your own override !")

    def _build_predefined(self, lines, network, station, origin, magnitude, NS, dataWindowRange):
        
        if len(NS) == 0:
            raise BadParameter()
        
        ta = min([ NS[k][-1] for k in NS ])
        tb = max([ NS[k][-1] for k in NS ])
        
        # We find some good defaults to reuse the method
        allowedGainCodes = list(set([ NS[k][2] + "." + NS[k][3][-1] for k in NS ]))
        
        targetSamplingRate = 50
        
        sps = {
            'L' : 1,
            'B' : 20,
            'H' : 100,
            'S' : 80,
            'E' : 500,
            'M' : 10,
            'C' : 1000
        }
        
        try:
            targetSamplingRate = max(list(set([ sps[NS[k][3][0]] for k in NS ])))
        except:
            print("\nFailed to callibrate sampling rate", file = sys.stderr)
        
        # This method is overriden in each implementation
        (z,n,e) = self._getChannelList(station, ta, targetSamplingRate, allowedGainCodes)
        
        EI = self.__build_event_dictionary(origin.time, origin.latitude, origin.longitude, origin.depth, magnitude.mag)
        SI = self.__build_station_dictionary(network.code, station.code, station.latitude, station.longitude, station.elevation)

        PI = None
        for k in [ 'P', 'Pn', 'pP', 'Pg', 'S', 'Sg', 'Sn' ]:
            if k not in NS: continue
            N,S,L,C,distance,azimuth,phasename,t = NS[k]
            PI = self.__build_pick_dictionary(phasename, t, 0.0, PI)

        lines.append((ta + dataWindowRange.min(),
                      tb + dataWindowRange.max(),
                      network.code,
                      station.code,
                      [z,n,e],
                      SI,
                      EI,
                      PI)
                     )

    def _build(self, lines, network, station, origin, magnitude, phaselist, phasename, targetSamplingRate, allowedGainCodes, timeRange):
        delta = geodetics.locations2degrees(station.latitude,
                                       station.longitude,
                                       origin.latitude,
                                       origin.longitude)

        (ta, slowness) = self.__find_arrivalTime(origin.time, delta, origin.depth, phaselist)

        # This method is overriden in each implementation
        (z,n,e) = self._getChannelList(station, ta, targetSamplingRate, allowedGainCodes)

        EI = self.__build_event_dictionary(origin.time, origin.latitude, origin.longitude, origin.depth, magnitude.mag)
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

    def _getOrigin(self, event):
        origin = event.preferred_origin()
        magnitude = event.preferred_magnitude()

        if origin is None:
            if len(event.origins) == 0:
                raise NextItem("Bad origin for event.")
            origin = event.origins[0]


        if magnitude is None:
            if len(event.magnitudes) == 0:
                raise NextItem("Bad Magnitude for event.")
            magnitude = event.magnitudes[0]

        if not isinstance(origin, ObsEvent.Origin): raise NextItem("Wrong origin for event.")
        if not isinstance(magnitude, ObsEvent.Magnitude): raise NextItem("Wrong magnitude for event.")

        return (origin, magnitude)

    def _resolve_phasenames(self, value):
        plist = []

        '''
        Define more groups as needed:
        '''
        groups = {
          "pgroup": [ "ttp" ], 
          "pall": [ "ttp" ]
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
        if 'Ot' not in phase:
            arrivals = self.__tworker.get_travel_times(depth / 1000.0, delta, phase_list=phase )
        else: # Handle the Ot flag to use Event Origin Time
            arrivals = [ AttribDict({
                'time' : 0.0,
                'ray_param_sec_degree' : 0.0
            })]

        if len(arrivals) == 0:
            raise NextItem("No phases selected for Depth: %s Distance: %s Phases: %s" % (delta, depth, phase))

        time = arrivals[0].time
        slowness = arrivals[0].ray_param_sec_degree

        return (t0 + time, slowness)

    def __build_event_dictionary(self, t0, originLatitude, originLongitude, originDepth, eventMagnitude):
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
                     ## Depth is to be in KM !
                     'depth': originDepth / 1000.0,
                     'magnitude': eventMagnitude
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
                       ## Elevation is to be in meters !
                       'elevation': stationElevation
                   }
        return AttribDict(stationinfo)

    def __build_pick_dictionary(self, phase, time, slowness, main = None):
        pickinfo = AttribDict({
                    'phase'    : phase,
                    'time'     : time,
                    'slowness' : slowness,
                    'others'   : None
                    })
        
        if main is not None:
            main['others'].append(pickinfo)
            return main
        
        pickinfo['others'] = []
        
        return pickinfo

    def _organize_by_station(self, lines):
        request = {}

        request["STATUS"] = Status()

        for line in lines:
            key = "%s.%s" % (line[2],line[3])
            try:
                clist = request[key]
            except KeyError:
                request[key] = []
                clist = request[key]
            clist.append(line)

        return request

    def _organize_by_event(self, lines):
        request= {}

        request["STATUS"] = Status()

        for line in lines:
            key = "%s" % (line[6]['eventId'])
            try:
                clist = request[key]
            except KeyError:
                request[key] = []
                clist = request[key]
            clist.append(line)

        return request

    def _fill_kwargsstation(self, t0, t1,
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

    def _fill_kwargsevent(self, t0, t1,
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

    def _check_param(self, t0, t1, targetSamplingRate, allowedGainCodes, timeRange, phasesOrPhaseGroup,
                        networkStationCodes, stationRestrictionArea,
                        eventRestrictionArea,
                        magnitudeRange,
                        depthRange,
                        distanceRange):

        if t0 is None: raise Exception("T0 should not be None")
        if t1 is None: raise Exception("T1 should not be None")

        if targetSamplingRate is None: raise Exception("targetSamplingRate should not be None")
        try:
            targetSamplingRate = float(targetSamplingRate)
        except ValueError:
            raise Exception("Invalid value for targetSamplingRate - expected float.")

        if allowedGainCodes is None or not isinstance(allowedGainCodes, list): raise Exception("allowedGainCodes should not be None - expected [ ... ]")
        if phasesOrPhaseGroup is None: raise Exception("phasesOrPhaseGroup should not be None - expected a phase name or group")

        if timeRange is None or not isinstance(timeRange, Range):
            raise Exception("Invalid TimeRange - expected Range Instance")

        if stationRestrictionArea is not None and not isinstance(stationRestrictionArea, AreaRange):
            raise Exception("Invalid StationRestrictionArea - expected AreaRange Instance")
        if eventRestrictionArea is not None and not isinstance(eventRestrictionArea, AreaRange):
            raise Exception("Invalid eventRestrictionArea - expected AreaRange Instance")

        if magnitudeRange is not None and not isinstance(magnitudeRange, Range):
            raise Exception("Invalid magnitudeRange - expected Range Instance")
        if depthRange is not None and not isinstance(depthRange, Range):
            raise Exception("Invalid depthRange - expected Range Instance")
        if distanceRange is not None and not isinstance(distanceRange, Range):
            raise Exception("Invalid distanceRange - expected Range Instance")

        if isinstance(t0, str): t0 = UTCDateTime(t0)
        if isinstance(t1, str): t1 = UTCDateTime(t1)

        return (t0, t1, targetSamplingRate, allowedGainCodes, timeRange, phasesOrPhaseGroup,
                    networkStationCodes, stationRestrictionArea,
                    eventRestrictionArea, magnitudeRange, depthRange,
                    distanceRange)

    ''' Public exported methods '''

    def setShowEvents(self, true_false):
        self._plotevents = bool(true_false)

    ''' Static methods have underscore '''

    @staticmethod
    def __x_list(data, fields, formats, validfields, formatrule, separator, destination):
        for f in formats:
            if f in formatrule:
                formatrule[f] = formats[f]
            else:
                raise Exception("Cannot change format for field %f" % f)

        for f in fields:
            if f not in validfields:
                raise Exception("Invalid field %f" % f)

        i = 1
        data_list = [ ]
        for e in data:
            ev = data[e]
            first = True
            data_fields = []
            for f in fields:
                if not first: print(separator, end = "", file=destination)
                if f == "#":
                    data_fields.append(i)
                    print(formatrule[f] % i, end = "", file=destination)
                else:
                    data_fields.append(ev[f])
                    print(formatrule[f] % ev[f],end = "", file=destination)
                first = False
            data_list.append(data_fields)
            print(file=destination)
            i += 1
        return data_list

    @staticmethod
    def stev_list(request, fields= ["#", "stationId", "slongitude", "slatitude", "selevation", "eventId", "etime", "elongitude", "elatitude", "edepth", "emagnitude"], separator = "\t", formats = { }, destination = sys.stdout):
        validfields = ["#", "stationId", "slongitude", "slatitude", "selevation", "eventId", "etime", "elongitude", "elatitude", "edepth", "emagnitude"]
        formatrule     = {
           "#": "%04d",
           "stationId": "%s",
           "slongitude": "%+9.4f",
           "slatitude": "%+9.4f",
           "selevation": "%6.2f",
           "eventId": "%s",
           "etime": "%s",
           "elongitude": "%+9.4f",
           "elatitude": "%+9.4f",
           "edepth": "%6.1f",
           "emagnitude": "%3.1f"
        }

        evs  = { }

        for key in request:
            if key == "STATUS": continue
            for line in request[key]:
                (_, _, _, _, _, SI, EI, _) = line
                if EI.eventId not in evs:
                    d = AttribDict()
                    d.eventId    = EI.eventId
                    d.etime      = EI.time
                    d.elongitude = EI.longitude
                    d.elatitude  = EI.latitude
                    d.emagnitude = EI.magnitude
                    d.edepth     = EI.depth
                    d.stationId  = SI.stationId
                    d.slongitude = SI.longitude
                    d.slatitude  = SI.latitude
                    d.selevation = SI.elevation
                    evs[(EI.eventId, SI.stationId)] = d

        return BaseBuilder.__x_list(evs, fields, formats, validfields, formatrule, separator, destination)

    @staticmethod
    def station_list(request, fields = ["#", "stationId", "longitude", "latitude", "elevation"], separator = "\t", formats = { }, destination = sys.stdout):
        validfields = ["#", "stationId", "longitude", "latitude", "elevation"]
        formatrule     = {
           "#": "%04d",
           "stationId": "%s",
           "longitude": "%+9.4f",
           "latitude": "%+9.4f",
           "elevation": "%6.2f",
        }

        sts  = { }

        for key in request:
            if key == "STATUS": continue
            for line in request[key]:
                (_, _, _, _, _, SI, _, _) = line
                if SI.stationId not in sts:
                    sts[SI.stationId] = SI

        return BaseBuilder.__x_list(sts, fields, formats, validfields, formatrule, separator, destination)

    @staticmethod
    def event_list(request, fields = ["#", "eventId", "time", "longitude", "latitude", "depth", "magnitude"], separator = "\t", formats = { }, destination = sys.stdout):
        validfields = ["#", "eventId", "time", "longitude", "latitude", "depth", "magnitude"]
        formatrule     = {
           "#": "%04d",
           "eventId": "%s",
           "time": "%s",
           "longitude": "%+9.4f",
           "latitude": "%+9.4f",
           "depth": "%6.1f",
           "magnitude": "%3.1f"
        }

        evs  = { }

        for key in request:
            if key == "STATUS": continue
            for line in request[key]:
                (_, _, _, _, _, _, EI, _) = line
                if EI.eventId not in evs:
                    evs[EI.eventId] = EI

        return BaseBuilder.__x_list(evs, fields, formats, validfields, formatrule, separator, destination)

    @staticmethod
    def filter_netsta(request, items, mode = "N.S"):
        '''
        Filter some stations/networks from the request. Mode is N.S, N or S to filter net/sta, net or stations
        Items is a list of codes to remove, following mode format.
        '''
        for evk in request:
            if evk == 'STATUS': continue
            ev = request[evk]
            lines = []
            for line in ev:
                if mode == "N.S":
                    k = "%s.%s" % (line[2], line[3])
                elif mode == "N":
                    k = "%s" % (line[2])
                elif mode == "S":
                    k = "%s" % (line[3])
                else:
                    raise BadParameter("Mode need to be N.S, N or S")
                if k not in items: continue
                lines.append(line)
            request[evk] = lines

        todelete = []
        for evk in request:
            if evk == 'STATUS': continue
            if len(request[evk]) == 0: todelete.append(evk)
        for evk in todelete:
            del request[evk]
        return request

    @staticmethod
    def filter_timewindow(request, start = None, end = None):
        if start is None and end is None: raise BadParameter("Need start or end parameters")
        for evk in request:
            if evk == 'STATUS': continue
            ev = request[evk]
            lines = []
            for line in ev:
                if start is not None and line[6]['time'] < start: continue
                if end is not None and line[6]['time'] > end: continue
                lines.append(line)
            request[evk] = lines

        todelete = []
        for evk in request:
            if evk == 'STATUS': continue
            if len(request[evk]) == 0: todelete.append(evk)
        for evk in todelete:
            del request[evk]
        return request

    @staticmethod
    def filter_channels(request, allowedChannels = "Z"):
        for evk in request:
            if evk == "STATUS" : continue
            ev = request[evk]
            lines = []
            for line in ev:
                items = list(filter(lambda x: allowedChannels in x[1], line[4]))
                line = (line[0], line[1], line[2], line[3], items, line[5], line[6], line[7])
                lines.append(line)
            request[evk] = lines
        return request

    @staticmethod
    def filter_netStationEvent(request, existing, useNetwork = False):
        keys = request.keys()
        ndel = 0
        for evk in keys:
            if evk == 'STATUS': continue
            ev = request[evk]
            lines = []
            for line in ev:
                if useNetwork:
                    t = ("%s.%s" % (line[2],line[3]), line[6].eventId[:-2])
                else:
                    t = (line[3], line[6].eventId[:-2])
                if t in existing:
                    ndel += 1
                    continue
                lines.append(line)
            if lines:
                request[evk] = lines
            else:
                del request[evk]
        print("  %d lines were removed" % ndel)
        return request

    @staticmethod
    def request_stats(request):
        ng = 0
        nl = { }
        
        start = {}
        end = {}
                
        for key in request:
            if key == "STATUS": continue
            ng += 1
            start[key] = []
            end[key] = []
            nl[key] = 0
            for line in request[key]:
                start[key].append(line[0])
                end[key].append(line[1])
                nl[key] += 1

        print("Total of ng=%d" % ng)
        print("Total of nl=%d" % sum(nl.values()))
        for k in nl:
            print("%s from %s to %s" % (k, min(start[k]), max(end[k]) ))
        print("Min Overall Date is=%s" % min(min(start.values())))
        print("Max Overall Date is=%s" % max(max(start.values())))
        
    @staticmethod
    def show_request(request, compact = True):
        def fv(v):
            if isinstance(v, float): return ("%7.1f" if not compact else "%.1f") % v
            return v

        def showkvd(title, it):
            if compact:
                print("%s " % (title), end="")
            else:
                print(title)

            for k,v in it.items():
                if compact:
                    print("%s: %s, " % (k,fv(v)), end="")
                else:
                    print("   %-10s: %s" % (k,fv(v)))

            if compact: print("")

        if request == None:
            raise Exception("Request cannot be None")

        print("\nRequent Entity of %d Line%s\n" % ((len(request) - 1), "s" if len(request) > 2 else ""))

        if "STATUS" in request:
                request["STATUS"].show()
                print("")

        for key in request:
            if key == "STATUS": continue
            print(key)
            for line in request[key]:
                print(" Start - End: %s - %s\n Netowrk Station: %s %s\n Streams: %s" % (line[0], line[1], line[2], line[3], line[4]))
                showkvd("  Station Attributes:", line[5])
                showkvd("  Event Attributes:", line[6])
                showkvd("  Pick Attributes:", line[7])
                print("")
        return

    @staticmethod
    def map_request(request, add_lines = False, enhance = None, showonly = None):
        from matplotlib import pyplot as plt
        from mpl_toolkits.basemap import Basemap

        keys = [ "ALL" ]
        for k in request:
            if k == "STATUS": continue
            keys.append(k)

        m = Basemap(projection='robin', lon_0=0, resolution='c')

        m.fillcontinents(color='gray',lake_color='white')
        m.drawmapboundary(fill_color='white')
        m.drawmeridians(range(0, 360, 30))
        m.drawparallels(range(-90, 90, 30))

        evn, stn, evs, sts = [ ], [ ], [ ], [ ]
        for k in request:
            if k == "STATUS": continue
            for item in request[k]:
                st = item[5]
                ev = item[6]
                eid = str(ev['eventId'])
                sid = str(st['stationId'])
                if showonly is not None and ((eid != showonly) and (sid != showonly)): continue
                elon = float(ev['longitude'])
                elat = float(ev['latitude'])
                slon = float(st['longitude'])
                slat = float(st['latitude'])
                ex, ey = m(elon, elat)
                sx, sy = m(slon, slat)
                if eid not in evn:
                    evn.append(eid)
                    evs.append((ex,ey))
                if sid not in stn:
                    stn.append(sid)
                    sts.append((sx,sy))
                if add_lines:
                    m.plot([ex, sx], [ey, sy], linewidth=1, color='k')

        for ex,ey in evs:
            m.plot(ex, ey, 'b*', markersize = 16, color = 'g')
        
        for sx,sy in sts:
            m.plot(sx, sy, 'b^', markersize = 16, color = 'r')
        
        if enhance is not None:
            for item in [ enhance ] if isinstance(enhance, str) else enhance:
                if item in stn:
                    ex, ey = sts[stn.index(item)]
                    m.plot(ex, ey, '^', markersize = 16, color = 'k')
                if item in evn:
                    ex, ey = evs[evn.index(item)]
                    m.plot(ex, ey, '*', markersize = 16, color = 'k')
        
        plt.title('Request Container')
        plt.show()

    @staticmethod
    def load_request(filename):
        if filename is None or not os.path.isfile(filename):
            raise Exception("Cannot read file, %s" % filename)

        try:
            iofile = open(filename, "rb")
            request = pickle.load(iofile)
            iofile.close()
        except:
            request = None

        return request

    @staticmethod
    def save_request(filename, request, overwrite = False):
        if filename is None:
            raise Exception("Filename should not be empty")

        if os.path.isfile(filename) and not overwrite:
            raise Exception("Will not overwrite file, %s" % filename)

        if request is None:
            raise Exception("Request cannot be empty")

        iofile = open(filename, "wb")
        pickle.dump(request, iofile)
        iofile.close()


if __name__ == "__main__":
    if not os.path.isfile("test_request.rq"):
        print("Please run prepare_tests.py first!")
        sys.exit(1)
    
    print("** Range **")
    print("All Depths = ",Range.ALLDEPTHS())
    print("All Dists  = ",Range.ALLDISTS())
    print("All Mags   = ",Range.ALLMAGS())
    print("")
    R = Range(0,1)
    print(R)
    print(" 0 in R= ", R.good(0), "= True")
    print(" 1 in R= ", R.good(1), "= True")
    print(" 2 in R= ", R.good(2), "= False")
    print("** End Range **")

    print("\n")

    print("** AreaRange **")
    print("World = ",AreaRange.WORLD())
    print("")
    R = AreaRange(0,1,0,1)
    print(R)
    print(" 0.5/0.5 in R= ", R.good(0.5,0.5))
    print(" 0/0 in R= ", R.good(0,0), "= True")
    print(" 0/1 in R= ", R.good(0,1), "= True")
    print(" 1/0 in R= ", R.good(1,0), "= True")
    print(" 1/1 in R= ", R.good(1,1), "= True")
    print(" 0/2 in R= ", R.good(0,2), "= False")
    print(" 2/0 in R= ", R.good(2,0), "= False")
    print(" 2/2 in R= ", R.good(2,2), "= False")
    print(" 0/-2 in R= ", R.good(0,-2), "= False")
    print(" -2/0 in R= ", R.good(-2,0), "= False")
    print(" -2/-2 in R= ", R.good(-2,-2), "= False")
    print("** End Range **")

    print("\n")

    print("** BaseBuilders **")
    
    # Assembly one request
    
    # Load
    
    try:
        BaseBuilder.load_request("test_request")
        print("Failed to refuse to load data")
    except:
        pass

    try:
        rq = BaseBuilder.load_request("test_request.rq")
    except:
        print("Failed to load a request")
        sys.exit(1)
    
    try:
        BaseBuilder.save_request("test_request.rq", rq, False)
        print("Failed to refuse to rewrite a file")
    except:
        pass

    # Lists
    
    print("\n** St/Ev List")
    BaseBuilder.stev_list(rq)
    print("**")
    
    print("\n** St List")
    BaseBuilder.station_list(rq)
    print("**")
    
    print("\n** Ev List")
    BaseBuilder.event_list(rq)
    print("**")
    
    print("\n** Show Request")
    BaseBuilder.show_request(rq, True)
    BaseBuilder.show_request(rq, False)
    print("**")
    
    # ~ BaseBuilder.map_request(rq) ### NOT WORKING
    
    BaseBuilder.request_stats(rq)
    
    # Filters
    
    BaseBuilder.filter_channels(rq, "Z")
    BaseBuilder.show_request(rq, True)
    
    BaseBuilder.filter_netsta(rq, [ "BL.AQDB" ])
    BaseBuilder.show_request(rq, True)
    
    BaseBuilder.filter_timewindow(rq, UTCDateTime("2010-01-01"), UTCDateTime("2011-01-01"))
    BaseBuilder.show_request(rq, True)
    
    # ~ BaseBuilder.filter_netStationEvent()  ### NOT TESTED
    
    print("** End BaseBuilders **")
