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

from obspy.clients import fdsn
from obspy.core import AttribDict
from obspy.core.utcdatetime import UTCDateTime
from obspy import geodetics
from obspy import taup
from seiscomp.arclink.manager import ArclinkManager
from obspy.core import event as ObsEvent

import pickle
import os

import socket

STATUS = AttribDict({
                     "unset": "unset",
                     
                     "requested": "requested",
                     "downloaded": "downloaded",
                     "saved": "saved",
                     
                     "temporary_error": "temporary_error",
                     "permanent_error": "permanent_error",
                     })

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

def unWrapNSLC(objs, archive = None, onlyShared = False):
    # unwrap lists of lists into arrays of tuples
    clist = []
    for (code, spam) in objs.items():
        for (start, obj) in spam.items():
            try:
                if archive and getattr(obj, "archive") != archive:
                    continue
                if onlyShared and getattr(obj, "shared") == False:
                    continue
            except:
                pass
            clist.append((code, start, obj))
    return clist

class Range(object):
    def __init__(self, minvalue, maxvalue):
        self._min = min([minvalue,maxvalue]) if minvalue != None and maxvalue != None else minvalue
        self._max = max([minvalue,maxvalue]) if minvalue != None and maxvalue != None else maxvalue

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
    def ALLMAGS():
        return Range(0.0, 10.0)

    @staticmethod
    def ALLDEPTHS():
        return Range(0.0, 1000.0)

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

class BaseBuilder(object):
    def __init__(self):
        self._plotevents = False
        self.__tworker = taup.TauPyModel()

    '''
    Generic Methods that can be used in FDSN or ARCLINK modes
    '''

    def _getChannelList(self, station, t0, targetsps, instcode):
        raise Exception("Implement your own override !")

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
          "pgroup": [ "ttp" ]
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
        arrivals = self.__tworker.get_travel_times(depth / 1000.0, delta, phase_list=phase )

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

    def __build_pick_dictionary(self, phase, time, slowness):
        pickinfo = {
                    'phase': phase,
                    'time': time,
                    'slowness': slowness
                    }
        return AttribDict(pickinfo)

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
        except ValueError,e:
            raise Exception("Invalid value for targetSamplingRate - expected float:\n  %s." % e.message)

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
        for e in data:
            ev = data[e]
            first = True
            for f in fields:
                if not first: print(separator, end = "", file=destination)
                if f == "#":
                    print(formatrule[f] % i, end = "", file=destination)
                else:
                    print(formatrule[f] % ev[f],end = "", file=destination)
                first = False
            print(file=destination)
            i += 1

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
        
        BaseBuilder.__x_list(evs, fields, formats, validfields, formatrule, separator, destination)

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
        
        BaseBuilder.__x_list(sts, fields, formats, validfields, formatrule, separator, destination)

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
        
        BaseBuilder.__x_list(evs, fields, formats, validfields, formatrule, separator, destination)

    @staticmethod
    def filter_channels(request, allowedChannels = "Z"):
        for evk in request:
            ev = request[evk]
            lines = []
            for line in ev:
                items = filter(lambda x: allowedChannels in x[1], line[4])
                line = (line[0], line[1], line[2], line[3], items, line[5], line[6], line[7])
                lines.append(line)
            request[evk] = lines
        return request

    @staticmethod
    def filter_netStationEvent(request, existing):
        keys = request.keys()
        for evk in keys:
            ev = request[evk]
            lines = []
            for line in ev:
                t = ("%s.%s" % (line[2],line[3]), line[6].eventId)
                if t in existing: continue
                lines.append(line)
            if lines:
                request[evk] = lines
            else:
                del request[evk]
        return request

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

            for k,v in it.iteritems():
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
    def map_request(request, add_lines = False):
        from matplotlib import pyplot as plt
        from mpl_toolkits.basemap import Basemap

        keys = [ "ALL" ]
        for k in request:
            if k == "STATUS": continue
            keys.append(k)

        def plot(showkey):
            if showkey == "ALL":
                showkey = None

            m = Basemap(projection='robin', lon_0=0, resolution='c')
            
            m.fillcontinents(color='gray',lake_color='white')
            m.drawmapboundary(fill_color='white')
            m.drawmeridians(range(0, 360, 30))
            m.drawparallels(range(-90, 90, 30))
            
            evn, stn, evs, sts = [ ], [ ], [ ], [ ]
            for k in request:
                if k == "STATUS": continue
                if showkey is not None and k != showkey: continue
                for item in request[k]:
                    st = item[5]
                    ev = item[6]
                    eid = ev['eventId']
                    sid = st['stationId']
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

            plt.title('Request Container')
            plt.show()

        try:
            from ipywidgets import interact
            interact(plot, showkey = keys);
        except ImportError:
            plot( "ALL" )

    @staticmethod
    def load_request(filename):
        if filename is None or not os.path.isfile(filename):
            raise Exception("Cannot read file, %s" % filename)

        try:
            iofile = file(filename, "r")
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

        iofile = file(filename, "w")
        pickle.dump(request, iofile)
        iofile.close()

class ArcLinkFDSNBuilder(BaseBuilder):
    '''Build a request using a FDSN server for events and an ArcLink for stations metadata'''
    def __init__(self, fdsn_event_url, arclink_url):
        BaseBuilder.__init__(self)
        (host,port,user) = arclink_url.strip().replace("/",":").split(":")
        self.__arclink_manager = ArclinkManager("%s:%s" % (host,port), user)
        self.__fdsn_client = fdsn.Client(fdsn_event_url)

    '''
    ArcLink specific Methods
    '''

    def __choose(self, item, location, channel, targetsps, instcode):

        # Build a newitem for the current channel
        sps = channel.sampleRateNumerator / channel.sampleRateDenominator
        newitem = (location.code, channel.code, sps, channel.azimuth, channel.dip, abs(sps - targetsps))

        # Check that the instrument code is allowed to select
        if channel.code[1] not in instcode: return item

        # If we did not select yet return the current
        if item[0] == None: return newitem

        # Decided based on SPS first
        if newitem[5] == item[5]:
            # And based on the Channel Instrument Code
            if instcode.index(newitem[1][1]) < instcode.index(item[1][1]):
                return newitem
        elif newitem[5] < item[5]:
            return newitem

        return item

    def _getChannelList(self, station, t0, targetsps, instcode):

        ## Choose will build something like: (loca, chan, sps , az  , dip , dsps)
        z = (None, None, None, None, None, None)
        n = (None, None, None, None, None, None)
        e = (None, None, None, None, None, None)

        for (_, _, location) in unWrapNSLC(station.sensorLocation):
            for (_, _, channel) in unWrapNSLC(location.stream):
                if t0 < channel.start or (channel.end is not None and t0 > channel.end):
                    continue
    
                if channel.code[-1] == "Z":
                    z = self.__choose(z, location, channel, targetsps, instcode)
                elif channel.code[-1] == "N" or channel.code[-1] == "1":
                    n = self.__choose(n, location, channel, targetsps, instcode)
                elif channel.code[-1] == "E" or channel.code[-1] == "2":
                    e = self.__choose(e, location, channel, targetsps, instcode)
    #             else:
    #                 print("Unknow channel %s.%s" % (channel.location_code, channel.code), file=sys.stderr)

        if z[1] is None or n[1] is None or e[1] is None:
            raise NextItem("No Z, N or E channel(s) found")

        # Return only (location, channel, azimuth, dip)
        z = (z[0], z[1], z[3], z[4])
        n = (n[0], n[1], n[3], n[4])
        e = (e[0], e[1], e[3], e[4])

        return (z,n,e)

    '''
    Request Builder Methods
    '''

    def eventBased(self, t0, t1, targetSamplingRate, allowedGainList, dataWindowRange, phasesOrPhaseGroupList,
                   eventRestrictionArea,
                   magnitudeRange,
                   depthRange,

                   networkStationList = None,
                   stationRestrictionArea = None,
                   distanceRange = None
                   ):
        '''Search initially by events, later, find stations related to the events time and distance specified.'''

        # List of request lines
        lines = []

        (phasename, phaselist) = self._resolve_phasenames(phasesOrPhaseGroupList)
        print("Searching using: %s %s" % (phasename, phaselist), file=sys.stderr)

        if networkStationList is None:
            networkStationList = [ "*.*" ]

        (t0, t1, targetSamplingRate, allowedGainList, dataWindowRange, phasesOrPhaseGroupList, 
        networkStationList, stationRestrictionArea,
        eventRestrictionArea, magnitudeRange, depthRange,
        distanceRange) = self._check_param(t0, t1, targetSamplingRate, allowedGainList, dataWindowRange, phasesOrPhaseGroupList, 
                                          networkStationList, stationRestrictionArea,
                                          eventRestrictionArea,
                                          magnitudeRange,
                                          depthRange,
                                          distanceRange)

        kwargsevent = self._fill_kwargsevent(t0, t1,
                                              eventRestrictionArea,
                                              magnitudeRange,
                                              depthRange,
                                              None, None, None)

        try:
            events = self.__fdsn_client.get_events(**kwargsevent)
            if self._plotevents:
                events.plot()
            print("Found %d events." % len(events), file=sys.stderr)
        except fdsn.header.FDSNException:
            print("No events found for the given parameters.", file=sys.stderr)
            return None

        # Event loop
        for event in events:
            try:
                (origin, magnitude) = self._getOrigin(event)
            except NextItem,e:
                print("Skipping Origin: %s" % str(e), file=sys.stderr)
                continue

            print("Working on origin: %s" % str(origin.time), file=sys.stderr)

            for code in networkStationList:
                (net, sta) = code.split(".")
                inventory = self.__arclink_manager.get_inventory(net, sta, "*", "*", (origin.time - 86400).datetime, (origin.time + 86400).datetime)

                if inventory is None or len(inventory.network) == 0:
                    print("No stations for pattern: %s.%s" % (net, sta), file=sys.stderr)
                    continue

                print("\n Stations for pattern: %s.%s" % (net,sta), file=sys.stderr)
                # Station loop
                for (_, _, network) in  unWrapNSLC(inventory.network):
                    for (_, _, station) in  unWrapNSLC(network.station):
                        # ! MISSING CHECKS
                        try:
                            print("  Working on station %s.%s " % (network.code, station.code), end="", file=sys.stderr)
                            self._build(lines,
                                         network, station, origin, magnitude,
                                         phaselist, phasename,
                                         targetSamplingRate, allowedGainList, dataWindowRange)
                            print("OK!", file=sys.stderr)
                        except NextItem, e:
                            print("\n  Skipping: %s" % str(e), file=sys.stderr)

        request = self._organize_by_event(lines)

        return request

    def stationBased(self, t0, t1, targetSamplingRate, allowedGainList, dataWindowRange, phasesOrPhaseGroupList,
                    networkStationList,
                    stationRestrictionArea,

                    eventRestrictionArea = None,
                    magnitudeRange = None,
                    depthRange = None,
                    distanceRange = None):
        '''Search initially by stations, later, find events related to the stations time and distance specified.'''

        # List of request lines
        lines = []

        (phasename, phaselist) = self._resolve_phasenames(phasesOrPhaseGroupList)
        print("Searching using: %s %s" % (phasename, phaselist), file=sys.stderr)

        (t0, t1, targetSamplingRate, allowedGainList, dataWindowRange, phasesOrPhaseGroupList, 
        networkStationList, stationRestrictionArea,
        eventRestrictionArea, magnitudeRange, depthRange,
        distanceRange) = self._check_param(t0, t1, targetSamplingRate, allowedGainList, dataWindowRange, phasesOrPhaseGroupList, 
                                          networkStationList, stationRestrictionArea,
                                          eventRestrictionArea,
                                          magnitudeRange,
                                          depthRange,
                                          distanceRange)

        # Check if network/station code is a string -> list
        if isinstance(networkStationList, str):
            networkStationList = [ networkStationList ]

        ## Start the Loop
        # On the given station patterns
        for code in networkStationList:
            (net, sta) = code.split(".")
            inventory = self.__arclink_manager.get_inventory(net, sta, "*", "*", t0.datetime, t1.datetime)

            if inventory is None or len(inventory.network) == 0:
                print("No stations for pattern: %s.%s" % (net, sta), file=sys.stderr)
                continue

            print("Stations for pattern: %s.%s" % (net,sta), file=sys.stderr)
            # Station Loop
            for (_, _, network) in  unWrapNSLC(inventory.network):
                for (_, _, station) in  unWrapNSLC(network.station):
                    print("\n Working on station %s.%s" % (network.code, station.code), file=sys.stderr)

                    # Start Build the parameters to pass on the Event Client
                    start = UTCDateTime(station.start) if station.start else t0
                    end = UTCDateTime(station.end) if station.end else t1
                    kwargsevent = self._fill_kwargsevent(max([start, t0]),
                                                        min([end, t1]),
                                                        eventRestrictionArea,
                                                        magnitudeRange,
                                                        depthRange,
                                                        station.latitude,
                                                        station.longitude,
                                                        distanceRange)

                    try:
                        events = self.__fdsn_client.get_events(**kwargsevent)
                    except fdsn.header.FDSNException,e:
                        events = None

                    if events is None or len(events) < 0:
                        print("  No Events Found.",file=sys.stderr)
                        continue

                    # Event loop
                    for event in events:
                        try:
                            (origin, magnitude) = self._getOrigin(event)
                            print("  Working on origin: %s" % str(origin.time), end="", file=sys.stderr)
                            self._build(lines,
                                       network, station, origin, magnitude,
                                       phaselist, phasename,
                                       targetSamplingRate, allowedGainList, dataWindowRange)
                            print("OK!", file=sys.stderr)
                        except NextItem,e:
                            print("  Skipping: %s" % str(e), file=sys.stderr)
                            continue

        request = self._organize_by_station(lines)

        return request

class FDSNBuilder(BaseBuilder):
    '''Build a request using a FDSN server for events and another for stations metadata'''
    def __init__(self, fdsn_event_url, fdsn_station_url = None):
        BaseBuilder.__init__(self)

        d = False

        if isinstance(fdsn_event_url, str):
            self.e_fdsn_client = fdsn.Client(fdsn_event_url, debug = d)
        elif isinstance(fdsn_event_url, fdsn.client):
            self.e_fdsn_client = fdsn_event_url
        else:
            raise BadParameter("Invalid event_serverorurl object, expected String address of fdsnClient Class")


        if fdsn_station_url == None:
            print("Same Server for Station as Event")
            self.s_fdsn_client = self.e_fdsn_client
        else:
            if isinstance(fdsn_station_url, str):
                self.s_fdsn_client = fdsn.Client(fdsn_station_url, debug = d)
            elif isinstance(fdsn_station_url, fdsn.client):
                self.s_fdsn_client = fdsn_station_url
            else:
                raise BadParameter("Invalid station_serverorurl object, expected String address of fdsnClient Class")

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
        elif newitem[5] < item[5]:
            return newitem

        return item

    def _getChannelList(self, station, t0, targetsps, instcode):

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
#                 print("Unknow channel %s.%s" % (channel.location_code, channel.code), file=sys.stderr)

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

    def stationBased(self, t0, t1, targetSamplingRate, allowedGainList, dataWindowRange, phasesOrPhaseGroupList,
                    networkStationList,
                    stationRestrictionArea,

                    eventRestrictionArea = None,
                    magnitudeRange = None,
                    depthRange = None,
                    distanceRange = None):
        '''Search initially by stations, later, find events related to the stations time and distance specified.'''

        (t0, t1, targetSamplingRate, allowedGainList, dataWindowRange, phasesOrPhaseGroupList, 
         networkStationList, stationRestrictionArea,
         eventRestrictionArea, magnitudeRange, depthRange,
         distanceRange) = self._check_param(t0, t1, targetSamplingRate, allowedGainList, dataWindowRange, phasesOrPhaseGroupList, 
                                          networkStationList, stationRestrictionArea,
                                          eventRestrictionArea,
                                          magnitudeRange,
                                          depthRange,
                                          distanceRange)

        (phasename, phaselist) = self._resolve_phasenames(phasesOrPhaseGroupList)
        print("Searching using: %s %s" % (phasename, phaselist), file=sys.stderr)

        
        # List of request lines
        lines = []

        # Start Build the parameters to pass on to the Client
        kwargsstation = self._fill_kwargsstation(t0,
                                                t1,
                                                stationRestrictionArea,
                                                None, None, None)

        # Check if network/station code is a string -> list
        if isinstance(networkStationList, str):
            networkStationList = [ networkStationList ]

        ## Start the Loop
        # On the given station patterns
        for code in networkStationList:
            (net, sta) = code.split(".")
            kwargsstation['net'] = net
            kwargsstation['sta'] = sta

            try:
                inventory = self.s_fdsn_client.get_stations(**kwargsstation)
            except:
                inventory = None

            if inventory is None or len(inventory.networks) == 0:
                print("No stations for pattern: %s.%s" % (net, sta), file=sys.stderr)
                continue

            # On networks found
            print("Stations for pattern: %s.%s" % (net,sta), file=sys.stderr)
            for network in inventory.networks:
                # On Stations found
                for station in network.stations:
                    print("\n Working on station %s.%s" % (network.code, station.code), file=sys.stderr)

                    evsdate = t0
                    if station.start_date and station.start_date > evsdate: evsdate = station.start_date

                    evedate = t1
                    if station.end_date and station.end_date < evedate: evedate = station.end_date

                    # Start Build the parameters to pass on the Event Client
                    kwargsevent = self._fill_kwargsevent(evsdate,
                                                          evedate,
                                                          eventRestrictionArea,
                                                          magnitudeRange,
                                                          depthRange,
                                                          station.latitude,
                                                          station.longitude,
                                                          distanceRange)

                    try:
                        events = "INVALID"
                        tryid = 1
                        while tryid < 3 and events == "INVALID":
                            try:
                                events = self.e_fdsn_client.get_events(**kwargsevent)
                            except socket.timeout:
                                print('Failed . try %d' % tryid)
                                tryid = tryid + 1
                                pass
                        if events == "INVALID": events = None
                    except fdsn.header.FDSNException,e:
                        events = None

                    # Event loop
                    if events is None or len(events) < 0:
                        print("  No Events Found.", file=sys.stderr)
                        continue

                    for event in events:
                        try:
                            (origin, magnitude) = self._getOrigin(event)

                            print("  Working on origin: %s" % str(origin.time), end="", file=sys.stderr)
                            self._build(lines,
                                       network, station, origin, magnitude,
                                       phaselist, phasename,
                                       targetSamplingRate, allowedGainList, dataWindowRange)
                            print("OK!", file=sys.stderr)
                        except NextItem,e:
                            print("  Skipping: %s" % str(e), file=sys.stderr)
                            continue

        request = self._organize_by_station(lines)

        return request

    def eventBased(self, t0, t1, targetSamplingRate, allowedGainList, dataWindowRange, phasesOrPhaseGroupList,
                   eventRestrictionArea,
                   magnitudeRange,
                   depthRange,

                   networkStationList = None,
                   stationRestrictionArea = None,
                   distanceRange = None
                   ):
        '''Search initially by events, later, find stations related to the events time and distance specified.'''
        lines = []

        (phasename, phaselist) = self._resolve_phasenames(phasesOrPhaseGroupList)
        print("Searching using: %s %s" % (phasename, phaselist), file=sys.stderr)

        if networkStationList is None:
            networkStationList = [ "*.*" ]

        (t0, t1, targetSamplingRate, allowedGainList, dataWindowRange, phasesOrPhaseGroupList, 
        networkStationList, stationRestrictionArea,
        eventRestrictionArea, magnitudeRange, depthRange,
        distanceRange) = self._check_param(t0, t1, targetSamplingRate, allowedGainList, dataWindowRange, phasesOrPhaseGroupList, 
                                          networkStationList, stationRestrictionArea,
                                          eventRestrictionArea,
                                          magnitudeRange,
                                          depthRange,
                                          distanceRange)

        kwargsevent = self._fill_kwargsevent(t0, t1,
                                              eventRestrictionArea,
                                              magnitudeRange,
                                              depthRange,
                                              None, None, None)

        try:
            events = self.e_fdsn_client.get_events(**kwargsevent)
            if self._plotevents:
                events.plot()
            print("Found %d events." % len(events), file=sys.stderr)
        except fdsn.header.FDSNException:
            print("No events found for the given parameters.", file=sys.stderr)
            return None

        # Event loop
        for event in events:
            try:
                (origin, magnitude) = self._getOrigin(event)
            except NextItem,e:
                print("Skipping Origin: %s" % str(e), file=sys.stderr)
                continue

            print("Working on origin: %s" % str(origin.time), file=sys.stderr)

            kwargsstation = self._fill_kwargsstation((origin.time - 86400),
                                                      (origin.time + 86400),
                                                      stationRestrictionArea,
                                                      origin.latitude,
                                                      origin.longitude,
                                                      distanceRange)

            for code in networkStationList:
                (net, sta) = code.split(".")
                kwargsstation['net'] = net
                kwargsstation['sta'] = sta

                try:
                    inventory = self.s_fdsn_client.get_stations(**kwargsstation)
                except fdsn.header.FDSNException,e:
                    inventory = None

                if inventory is None or len(inventory.networks) == 0:
                    print("\n No stations for pattern: %s.%s" % (net, sta), file=sys.stderr)
                    continue

                print("\n Stations for pattern: %s.%s" % (net,sta), file=sys.stderr)
                # Event loop
                for network in inventory.networks:
                    for station in network.stations:
                        try:
                            print("  Working on station %s.%s " % (network.code, station.code), end="", file=sys.stderr)
                            self._build(lines,
                                       network, station, origin, magnitude,
                                       phaselist, phasename,
                                       targetSamplingRate, allowedGainList, dataWindowRange)
                            print("OK!", file=sys.stderr)
                        except NextItem, e:
                            print("\n  Skipping: %s" % str(e), file=sys.stderr)

        request = self._organize_by_event(lines)

        return request

if __name__ == "__main__":
#     rb = ArcLinkRequestBuilder("IRIS","seisrequest.iag.usp.br:18001:m.bianchi@iag.usp.br")
#     rb = FDSNBuilder("IRIS")

#     req = rb.load_request("request.pik")

#     Call the stationBased
#     req = rb.eventBased(t0 = UTCDateTime("2007-01-01"),
#                         t1 = UTCDateTime("2008-01-01"),
#                         targetSamplingRate = 20.0,
#                         allowedGainList = ["H", "L"],
#                         dataWindowRange = Range(-120, 600),
#                         phasesOrPhaseGroupList = "pgroup",
# 
#                         networkStationList = [ "TA.A*"],
#                         stationRestrictionArea = AreaRange(-150.0, -90.0, 15.0, 60.0),
# 
#                         eventRestrictionArea = AreaRange(-35.0, -10.0, -55.0, -60.0),
#                         magnitudeRange = Range(6.2, 9.0),
#                         depthRange = Range(0.0, 400.0),
#                         distanceRange = None
#                         )
# 
#     rb.show_request(req)

#     rb.save_request("request.pik", req)
#     rb.show_request(req)

#     rq = rb.load_request("xc-phase2")
#     rb.event_list(rq, separator="\t")
#     rb.station_list(rq, separator="\t")
        
        pass
