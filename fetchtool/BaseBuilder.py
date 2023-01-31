'''BaseBuilder Module for the fetchtool package.

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
''' Basic status information dictionary with conversion strings '''

MINUTE = 60
''' Minute constant in seconds'''

HOUR   = 60 * MINUTE
''' Hour constant in seconds'''

DAY    = 24 * HOUR
''' One day constant in seconds'''

class Status(object):
    '''A not so usefull status representation of the request.
    '''
    def __init__(self):
        self.level = STATUS.unset
        '''The current status indication'''
        self.comment = None
        '''Any added comment to the status'''

    def __iter__(self):
        return self

    def __len__(self):
        return 0

    def next(self):
        '''Used for iteration'''
        raise StopIteration

    def show(self):
        '''Print a status message'''
        print("Status Level: ",self.level, "[",self.comment,"]")


class BadParameter(Exception):
    ''' Exception to indicate that a wrong parameter was supplied by the user '''
    pass


class NextItem(Exception):
    ''' Exception to indicate that we fail find some information and the iteration must proceed to the next item in list '''
    pass


class Range(object):
    '''A generic interval (start/end) representation
    
    Parameter
    ---------
    minvalue : number
        The start value, or minimum value of the interval.
    maxvalue : number
        The end value, or minimum value of the interval.
    '''
    def __init__(self, minvalue, maxvalue):
        self._min = minvalue # min([minvalue,maxvalue]) if minvalue != None and maxvalue != None else minvalue
        self._max = maxvalue # max([minvalue,maxvalue]) if minvalue != None and maxvalue != None else maxvalue

    def min(self):
        '''The left value of Range
        
        Return
        ------
        number
            The left (or minimum) value of the Range
        '''
        return self._min

    def max(self):
        '''The right value of Range
        
        Return
        ------
        number
            The right (or maximum) value of the Range
        '''
        return self._max

    def good(self, value):
        '''Check that value is a good value. I.e. is inside the Range
        
        Parameters
        ----------
        value : number
            A value to check
        
        Return
        ------
        boolean
            True of False for the test
        '''
        
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
        '''All possibles distances on Earth
        
        Return
        ------
        AreaRange
            An AreaRange representing all possible distances (0,180) degrees
        '''
        return Range(0.0, 180.0)
    
    @staticmethod
    def ALLMAGS():
        '''All possibles magnitudes
        
        Return
        ------
        AreaRange
            An AreaRange representing all possible magnitudes (0,10)
        '''
        return Range(0.0, 10.0)

    @staticmethod
    def ALLDEPTHS():
        '''All possibles depths
        
        Return
        ------
        AreaRange
            An AreaRange representing all possible depths (0,1000) in km
        '''
        return Range(0.0, 1000.0)

    def __str__(self):
        return "Range from: {} to: {}".format(self.min(), self.max())


class AreaRange(object):
    '''A generic region (squared) representation
    
    Parameter
    ---------
    xmin : number
        Left (west) limit of the region.
    xmax : number
        Right (east) limit of the region.
    ymin : number
        Bottom (south) limit of the region.
    ymax : number
        Top (north) limit of the region.
    '''
    
    def __init__(self, xmin, xmax, ymin, ymax):
        self.x = Range(xmin, xmax)
        '''Internal Range to represent the x-range of the AreaRange'''
        self.y = Range(ymin, ymax)
        '''Internal Range to represent the y-range of the AreaRange'''

    def xmin(self):
        '''The left value of AreaRange
        
        Return
        ------
        number
            The left value of the AreaRange
        '''
        return self.x.min()

    def xmax(self):
        '''The right value of AreaRange
        
        Return
        ------
        number
            The right value of the AreaRange
        '''
        return self.x.max()

    def ymin(self):
        '''The bottom value of AreaRange
        
        Return
        ------
        number
            The bottom value of the AreaRange
        '''
        return self.y.min()

    def ymax(self):
        '''The top value of AreaRange
        
        Return
        ------
        number
            The top value of the AreaRange
        '''
        return self.y.max()

    def good(self, x, y):
        '''Check that x,y point is a good point. I.e. is inside the AreaRange
        
        Parameters
        ----------
        x : number
            The x-coordinate of the point to check
        y : number
            The y-coordinate of the point to check
        
        Return
        ------
        boolean
            True of False for the test
        '''
        return self.x.good(x) and self.y.good(y)

    @staticmethod
    def WORLD():
        '''AreaRange representing the world region
        
        Return
        ------
        AreaRange
            An AreaRange object representing the world
        '''
        return AreaRange(-180., 180., -90., 90)

    @staticmethod
    def BRAZIL():
        '''AreaRange representing the Brazil region
        
        Return
        ------
        AreaRange
            An AreaRange object representing Brazil
        '''
        return AreaRange(-75., -30., -35., 8.)

    def __str__(self):
        return 'AreaRange from {} to {} and {} to {}'.format(self.xmin(), self.xmax(), self.ymin(), self.ymax())


class BaseBuilder(object):
    '''This is the Builders Super Class
    
    All Builders should implement the BaseBuilder class
    '''
    def __init__(self):
        self._plotevents = False
        '''A global flag used by eventBased and stationBased to show or not events.'''
        
        self.__include_restricted = True
        '''A global flag used by self._fill_kwargsstation to include or not restricted inventory.'''
        
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

    def _build_raw_lines(self, lines, ta, tb, network, station, targetSamplingRate, allowedGainCodes):
        # This method is overriden in each implementation
        (z,n,e) = self._getChannelList(station, ta, targetSamplingRate, allowedGainCodes)

        SI = self.__build_station_dictionary(network.code, station.code, station.latitude, station.longitude, station.elevation)

        lines.append((ta,
                      tb,
                      network.code,
                      station.code,
                      [z,n,e],
                      SI,
                      None,
                      None)
                     )

        return

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
                  "level": "channel",
                  "includerestricted" : self.__include_restricted
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
        '''Show events obtained from server before continuing to build the request.
        
        Parameter
        ---------
        true_false : bool
            A True or False value that will activate or disable the Show Event feature.
        '''
        self._plotevents = bool(true_false)

    def setAllowRestricted(self, true_false):
        '''Enable of disable the use of restricted station while building the request
        
        Parameters
        ----------
        true_false : bool
            A flag to enable (true) disable (false) the use of restricted stations.
        '''
        self.__include_restricted = bool(true_false)

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
        '''Generate a printout list of events and stations in the request
        
        This method can be used to generate a print out list of the stations and events in the request.
        Its output can be customized using the fields parameters. Aditionaly it also
        returns a list with the output that can be consumed into your program.
        
        Parameters
        ----------
        request : request
            The input request
        fields: list
            A list of fields to be reported and returned.
        separator : str, default "\t"
            A character that will be used as field separator (like "," to make a csv)
        formats : dict, default {}
            A dictionary of formats to change each field number formating.
        destination : file, default sys.stdout
            An open file for destination of the print out.
        
        Return
        ------
        list
            A data list including the selected fields.
        '''
        
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
                if EI is None:
                    print("Request is for continous data - no event information is available", file = sys.stderr)
                    return []
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
        '''Generate a printout list of individual stations in the request
        
        This method can be used to generate a print out list of the stations in the request.
        Its output can be customized using the fields parameters. Aditionaly it also
        returns a list with the output that can be consumed into your program.
        
        Parameters
        ----------
        request : request
            The input request
        fields: list
            A list of fields to be reported and returned.
        separator : str, default "\t"
            A character that will be used as field separator (like "," to make a csv)
        formats : dict, default {}
            A dictionary of formats to change each field number formating.
        destination : file, default sys.stdout
            An open file for destination of the print out.
        
        Return
        ------
        list
            A data list including the selected fields.
        '''
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
        '''Generate a printout list of individual events in the request
        
        This method can be used to generate a print out list of the events in the request.
        Its output can be customized using the fields parameters. Aditionaly it also
        returns a list with the output that can be consumed into your program.
        
        Parameters
        ----------
        request : request
            The input request
        fields: list
            A list of fields to be reported and returned.
        separator : str, default "\t"
            A character that will be used as field separator (like "," to make a csv)
        formats : dict, default {}
            A dictionary of formats to change each field number formating.
        destination : file, default sys.stdout
            An open file for destination of the print out.
        
        Return
        ------
        list
            A data list including the selected fields.
        '''
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
                if EI is None:
                    print("Request is for continous data - no event information is available", file = sys.stderr)
                    return []
                if EI.eventId not in evs:
                    evs[EI.eventId] = EI

        return BaseBuilder.__x_list(evs, fields, formats, validfields, formatrule, separator, destination)

    @staticmethod
    def filter_netsta(request, items, mode = "N.S", operation = "keep"):
        '''Filter the request based on the items, mode and operation.
        
        This method can be used to filter IN (operation="keep") or filter OUT (operation="remove") request lines
        from request based on its network/station values. The mode parameter indicate if items is a list of networks,
        stations of network.station values. Items is always a list.
        
        Warnings
        --------
        Requests is filtered inplace.
        
        Parameters
        ----------
        request : request
            A request to process. This will be modified.
        items : list
            A list of items to keep in the request.
        mode : str
            A string to indicate the mode of operation, can be N.S, N or S to indicate Network and Stations, Network or Station.
        operation : str, default "keep"
            A parameter to indicate if we should keep the supplied values (operation="keep", default) or remove (operation="remove") from the request.
        
        Return
        ------
        request 
            The request
        '''
        
        if operation not in ["keep", "remove"]:
            raise BadParameter("Invalid operation mode, valid values are keep or remove")
        
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
                if operation == "keep" and k not in items: continue
                if operation == "remove" and k in items: continue
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
        '''Filter out request lines not overlaping with this datetime range
        
        This method will remove request lines from the request that do not overlap the
        indicated time interval. The time interval can be open to the left or right.
        
        Warnings
        --------
        Requests are filtered inplace.
        
        Parameters
        ----------
        request : request
            The input request to be filtered
        start : str or UTCDateTime, default None
            The start of the interval to consider. None value indicates that the interval is open to the left.
        end : str or UTCDateTime, default None
            The end of the interval to consider. None value indicates that the interval is open to the right.
        
        Return
        ------
        request
            The filtered request
        '''
        if start is None and end is None: raise BadParameter("Need start or end parameters")
        
        if start is not None:
            start = UTCDateTime(start)
        
        if end is not None:
            end = UTCDateTime(end)
        
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
        '''Filter channels to be requested.
        
        This method will filter OUT all channels that are not in the allowedChannels variable.
        
        Warnings
        --------
        Request is filtered inplace.
        
        Parameters
        ----------
        request : request
            The input request for filtering
        allowedChannels : str, default "Z"
            A string with all allowed channel codes. Examples are: "Z", "12", "NE"
        
        Return
        ------
        request
            The filtered request
        '''
        
        for evk in request:
            if evk == "STATUS" : continue
            ev = request[evk]
            lines = []
            for line in ev:
                items = list(filter(lambda x: x[1][-1] in allowedChannels, line[4]))
                if len(items) == 0:
                    continue
                line = (line[0], line[1], line[2], line[3], items, line[5], line[6], line[7])
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
    def filter_netStationEvent(request, existing, useNetwork = False, precision = 60.0):
        '''
        Filter out request lines that are already listed in the existing list.
        
        existing is a list of (Network.Station, EventId ) tuples. Network can 
        be ommited if useNetwork = False.
        
        eventIds encodes event origin times, when comparing a tolerance of 
        precision is assumed due to changes in event location between built
        requests.
        '''
        def __check__(one, existing):
            for node, oid in existing:
                if node == one[0] and abs(oid - one[1]) < precision: return True
            return False
        
        existing = [ (node, UTCDateTime.strptime(evid, "%y%m%d_%H%M%S")) for node, evid in existing ]
        
        keys = request.keys()
        ndel = 0
        for evk in keys:
            if evk == 'STATUS': continue
            ev = request[evk]
            lines = []
            for line in ev:
                if useNetwork:
                    t = ("%s.%s" % (line[2],line[3]), UTCDateTime.strptime(line[6].eventId, "%y%m%d_%H%M%S"))
                else:
                    t = (line[3], UTCDateTime.strptime(line[6].eventId, "%y%m%d_%H%M%S"))
                if __check__(t, existing):
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
        '''Show some stats from the request
        
        Parameter
        ---------
        request : request
            The input request
        
        '''
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
    def reqlen(request):
        '''Count the total number of lines in request
        
        Parameters
        ----------
        request : request
            The request
        
        Return
        ------
        int
            The number of lines of request in request
        '''
        return sum([ len(request[k]) for k in request if k != "STATUS" ])
    
    @staticmethod
    def show_request(request, compact = True):
        '''Write a summary of the  request to the screen
        
        Parameters
        ----------
        request : request
            A request
        compact : bool, default True
            Use a compact mode to save space on screen
        '''
        
        def fv(v):
            if isinstance(v, float): return ("%7.1f" if not compact else "%.1f") % v
            return v

        def showkvd(title, it):
            if compact:
                print("%s " % (title), end="")
            else:
                print(title)

            if it is not None:
                for k,v in it.items():
                    if compact:
                        print("%s: %s, " % (k,fv(v)), end="")
                    else:
                        print("   %-10s: %s" % (k,fv(v)))
            else:
                if compact:
                    print(" - n/a - ", end = "")
                else:
                    print(" - n/a - ")

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
    def map_request(request, add_lines = False, enhance = None, 
                    showonly = None, global_map = True, grids_on = True):
        '''Make a map using cartopy representing the request.
        
        It plots all individual stations and events in the request. Map can be customized from the input parameters.
        
        Parameter
        ---------
        request : request
            The input request to plot
        add_lines : bool, default False
            Show a line connecting events and stations on the request
        enhance : str, default None
            A list of stations ids and/or events ids to be enhanced on the map
        showonly : str, default None
            One eventid or one stationid to be shown.
        global_map : bool, default True
            Show a global map (True) or a map around stations and events (False)
        grids_on : bool, default True
            Enable or disable the display of grids on the map
        
        '''
        import warnings
        from cartopy import crs, feature as cfeature
        from matplotlib import pyplot as plt
        
        warnings.filterwarnings('ignore')

        keys = [ "ALL" ]
        for k in request:
            if k == "STATUS": continue
            keys.append(k)

        fig = plt.figure(figsize=(11, 8.5))
        
        proj = crs.Robinson(central_longitude=0.0)
        
        data_trans=crs.Geodetic()
        
        ax = plt.subplot(1, 1, 1, projection = proj)
        ax.coastlines()
        
        ax.add_feature(cfeature.BORDERS, linewidth=0.5, edgecolor='black')
        
        evn, stn, evs, sts = [ ], [ ], [ ], [ ]
        for k in request:
            if k == "STATUS": continue
            for item in request[k]:
                st = item[5]
                ev = item[6]
                
                eid = str(ev['eventId']) if ev is not None  else None
                sid = str(st['stationId'])

                if showonly is not None and ((eid != showonly) and (sid != showonly)): continue
                
                slon = float(st['longitude'])
                slat = float(st['latitude'])
                sx, sy = slon, slat #m(slon, slat)

                if ev is not None:
                    elon = float(ev['longitude'])
                    elat = float(ev['latitude'])
                    ex, ey = elon, elat #m(elon, elat)
                    
                    if eid not in evn:
                        evn.append(eid)
                        evs.append((ex,ey))
                    
                if sid not in stn:
                    stn.append(sid)
                    sts.append((sx,sy))

                if add_lines and ev is not None:
                    ax.plot([ex, sx], [ey, sy], linewidth=1, color='k', transform=data_trans)

        
        for ex,ey in evs:
            ax.plot(ex, ey, 'b*', markersize = 16, color = 'g', transform = data_trans)
        
        for sx,sy in sts:
            ax.plot(sx, sy, 'b^', markersize = 16, color = 'r', transform = data_trans)
        
        if enhance is not None:
            for item in [ enhance ] if isinstance(enhance, str) else enhance:
                if item in stn:
                    ex, ey = sts[stn.index(item)]
                    ax.plot(ex, ey, '^', markersize = 16, color = 'k', transform = data_trans)
                if item in evn:
                    ex, ey = evs[evn.index(item)]
                    ax.plot(ex, ey, '*', markersize = 16, color = 'k', transform = data_trans)
        
        if global_map:
            ax.set_global()
        
        if grids_on:
            ax.gridlines(draw_labels=True, dms=True, x_inline=False, y_inline=False)

        plt.title('Request Container')
        plt.show()

    @staticmethod
    def load_request(filename):
        '''Load a request from file
        
        Parameters
        ----------
        filename : str
            A filename to read request from
        Return
        ------
        request
            The request object loaded from filename
        '''
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
        '''Save a built request to file
        
        Parameters
        ----------
        filename : str
            The filename to write to.
        request : request
            The request to be written.
        overwrite : bool, default True
            If the could should overwrite the file if it exists.
        '''
        
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
