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
from BaseBuilder import BaseBuilder, NextItem, BadParameter
import sys,socket

'''
  Clients
'''
from obspy.clients import fdsn

'''
  Tools
'''
from obspy import UTCDateTime
from obspy.core import AttribDict
from obspy import geodetics

'''
  Builders Definition
'''
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
    @staticmethod
    def _choose(item, channel, targetsps, instcode):
        def which(item, instcode):
            loc     = item[0]
            channel = item[1]
            for i,code in enumerate(instcode):
                codel,codec = ("*",code)
                if "." in code: codel,codec = code.split(".")
                
                if codel == "*":
                     if codec == channel[1]:
                         return i
                else:
                    if codel == loc and channel[1] == codec: return i
            
            return -1

        # Build a newitem for the current channel
        newitem = (channel.location_code, channel.code, channel.sample_rate, channel.azimuth, channel.dip, abs(channel.sample_rate - targetsps))

        # Check that the instrument code is allowed to select
        if which(newitem, instcode) == -1: return item

        # If we did not select yet return the current
        if item[0] == None: return newitem

        # Decided based on SPS first
        if newitem[5] == item[5]:
            # And based on the Channel Instrument Code
            if which(newitem, instcode) < which(item, instcode):
                return newitem
        elif newitem[5] < item[5]:
            return newitem

        return item

    @staticmethod
    def _getChannelList(station, t0, targetsps, instcode):

        ## Choose will build something like: (loca, chan, sps , az  , dip , dsps)
        z = (None, None, None, None, None, None)
        n = (None, None, None, None, None, None)
        e = (None, None, None, None, None, None)

        for channel in station.channels:
            if t0 < channel.start_date or (channel.end_date is not None and t0 > channel.end_date):
                continue

            if channel.code[-1] == "Z":
                z = FDSNBuilder._choose(z, channel, targetsps, instcode)
            elif channel.code[-1] == "N" or channel.code[-1] == "1":
                n = FDSNBuilder._choose(n, channel, targetsps, instcode)
            elif channel.code[-1] == "E" or channel.code[-1] == "2":
                e = FDSNBuilder._choose(e, channel, targetsps, instcode)
#             else:
#                 print("Unknow channel %s.%s" % (channel.location_code, channel.code), file=sys.stderr)

        if z[0] is None or n[0] is None or e[0] is None:
            raise NextItem("No Z, N or E channel(s) found considering loc/gain=%s" % (instcode))

        if z[0] != n[0] != e[0]:
            raise NextItem("Wow, algorith selected different location codes for different orientation!")

        if z[2] != n[2] != e[2]:
            raise NextItem("Wow, algorith selected different sampling codes for different orientation!")

        # Return only (location, channel, azimuth, dip)
        z = (z[0], z[1], z[3], z[4])
        n = (n[0], n[1], n[3], n[4])
        e = (e[0], e[1], e[3], e[4])

        return (z,n,e)

    '''
    Request Builder Methods
    '''

    def stationBased(self, t0, t1, targetSamplingRate, allowedLocGainList, dataWindowRange, phasesOrPhaseGroupList,
                    networkStationList,
                    stationRestrictionArea,

                    eventRestrictionArea = None,
                    magnitudeRange = None,
                    depthRange = None,
                    distanceRange = None):
        '''Search initially by stations, later, find events related to the stations time and distance specified.'''

        (t0, t1, targetSamplingRate, allowedLocGainList, dataWindowRange, phasesOrPhaseGroupList,
         networkStationList, stationRestrictionArea,
         eventRestrictionArea, magnitudeRange, depthRange,
         distanceRange) = self._check_param(t0, t1, targetSamplingRate, allowedLocGainList, dataWindowRange, phasesOrPhaseGroupList,
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
                        while tryid < 5 and events == "INVALID":
                            try:
                                print("Getting events ..")
                                events = self.e_fdsn_client.get_events(**kwargsevent)
                            except socket.timeout:
                                print('Failed . try %d' % tryid)
                                tryid = tryid + 1
                                pass
                            except socket.error:
                                print('Failed . try %d' % tryid)
                                tryid = tryid + 1
                                pass
                        if events == "INVALID": events = None
                    except fdsn.header.FDSNException as e:
                        print("FAILED:", e)
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
                                       targetSamplingRate, allowedLocGainList, dataWindowRange)
                            print("OK!", file=sys.stderr)
                        except NextItem as e:
                            print("  Skipping: %s" % str(e), file=sys.stderr)
                            continue

        request = self._organize_by_station(lines)

        return request

    def eventBased(self, t0, t1, targetSamplingRate, allowedLocGainList, dataWindowRange, phasesOrPhaseGroupList,
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

        (t0, t1, targetSamplingRate, allowedLocGainList, dataWindowRange, phasesOrPhaseGroupList,
        networkStationList, stationRestrictionArea,
        eventRestrictionArea, magnitudeRange, depthRange,
        distanceRange) = self._check_param(t0, t1, targetSamplingRate, allowedLocGainList, dataWindowRange, phasesOrPhaseGroupList,
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
            except NextItem as e:
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
                except fdsn.header.FDSNException as e:
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
                                       targetSamplingRate, allowedLocGainList, dataWindowRange)
                            print("OK!", file=sys.stderr)
                        except NextItem as e:
                            print("\n  Skipping: %s" % str(e), file=sys.stderr)

        request = self._organize_by_event(lines)

        return request

class CSVBuilder(BaseBuilder):
    '''
    Build a request from a Text File for Event and a FDSN/ArcLink server for stations.
    File format is:
    
    <Origin Time> <Longitude> <Latitude> <Depth> <Magnitude> [...]
    
    Separators and also be [ ], [,] or even [;].
    '''

    def __init__(self, time_lon_lat_dep_mag_file, arclink_or_fdsn, debug = False):
        BaseBuilder.__init__(self)
        
        self._evfile = time_lon_lat_dep_mag_file
        if arclink_or_fdsn.find("http://") == -1:
            (host,port,user) = arclink_or_fdsn.strip().replace("/",":").split(":")
            self._client = ArclinkManager("%s:%s" % (host,port), user)
            self._client_type = 'arclink'
        else:
            self._client = fdsn.Client(arclink_or_fdsn, debug = debug)
            self._client_type = 'fdsn'
        
        return

    '''
    CSV specific Methods
    '''
    def _getChannelList(self, station, t0, targetsps, instcode):
        
        if self._client_type == "arclink":
            return ArcLinkFDSNBuilder._getChannelList(station, t0, targetsps, instcode)
        
        if self._client_type == "fdsn":
            return FDSNBuilder._getChannelList(station, t0, targetsps, instcode)
        
        raise Exception("Bad Client!")
    
    def __loadcsv(self, t0, t1, evr, mr, dr, station = None, distanceRange = None):
        event     = AttribDict()
        magnitude = AttribDict()
        
        if station is not None:
            if self._client_type == "arclink":
                t0 = station.start
                t1 = station.end
            elif self._client_type == "fdsn":
                t0 = station.start_date
                t1 = station.end_date
        
        with open(self._evfile, "r") as fio:
            for line in fio:
                '''
                Parsers
                '''
                line = line.replace(";", " ")
                line = line.replace(",", " ")
                it = line.split()
                if len(it) < 5: 
                    raise Exception("Bad file %s given -- only %d columns found on '%s'" % (self._evfile, len(it), line))
                try:
                    lon = float(it[1]); lat = float(it[2]); dep = float(it[3]); mag = float(it[4]); time = UTCDateTime(it[0])
                except ValueError:
                    raise Exception("Bad file %s given -- only %d columns found on '%s -- no conversion was possible.'" % (self._evfile, len(it), line))
                
                '''
                Checking constraints
                '''
                if t0 is not None and time < t0: continue
                if t1 is not None and time > t1: continue
                if evr is not None and evr.good(lon, lat) == False: continue
                if mr is not None and mr.good(mag) == False: continue
                if dr is not None and dr.good(dep) == False: continue
                
                # Build the information into AttribDict
                ##
                event.time = time
                event.latitude = lat
                event.longitude = lon
                event.depth = dep
                magnitude.mag = mag
                
                '''
                Yield Data
                '''
                yield (event, magnitude)

        raise StopIteration()
    
    def __loadstArclink(self, net, sta, t0, t1, stationRestrictionArea, origin = None, distanceRange = None):
        if origin is not None:
            t0 = origin.time - 86400
            t1 = origin.time + 86400
        
        inventory = self._client.get_inventory(net, sta, "*", "*", t0.datetime, t1.datetime)
        
        if inventory is None or len(inventory.network) == 0:
            print("\n No stations for pattern: %s.%s" % (net, sta), file=sys.stderr)
            raise StopIteration("No Stations.")
        
        # Station loop
        for (_, _, network) in  unWrapNSLC(inventory.network):
            for (_, _, station) in  unWrapNSLC(network.station):
                '''
                Checks
                '''
                if origin is not None and distanceRange is not None:
                    delta = geodetics.locations2degrees(station.latitude,
                                               station.longitude,
                                               origin.latitude,
                                               origin.longitude)
                    if not distanceRange.good(delta): continue
                if stationRestrictionArea is not None:
                    if not stationRestrictionArea.good(station.longitude, station.latitude):
                        continue
                '''
                Yield Data
                '''
                yield network, station
        
        raise StopIteration()
    
    def __loadstFDSN(self, net, sta, t0, t1, stationRestrictionArea, origin = None, distanceRange = None):
        if origin is not None:
            t0 = origin.time - 86400
            t1 = origin.time + 86400
        
        kwargsstation = self._fill_kwargsstation(t0,
                                                t1,
                                                stationRestrictionArea,
                                                origin.latitude if origin is not None else None,
                                                origin.longitude if origin is not None else None,
                                                distanceRange)
        kwargsstation['net'] = net
        kwargsstation['sta'] = sta
        
        try:
            inventory = self._client.get_stations(**kwargsstation)
        except fdsn.header.FDSNException:
            inventory = None
        
        if inventory is None or len(inventory.networks) == 0:
            print("\n No stations for pattern: %s.%s" % (net, sta), file=sys.stderr)
            raise StopIteration("No Stations.")
        
        for network in inventory.networks:
            for station in network.stations:
                '''
                Yield Data
                '''
                yield network, station
        
        raise StopIteration()

    '''
    Request Builder Methods
    '''
    def eventBased(self, t0, t1, targetSamplingRate, allowedLocGainList, dataWindowRange, phasesOrPhaseGroupList,
                   eventRestrictionArea,
                   magnitudeRange,
                   depthRange,

                   networkStationList = None,
                   stationRestrictionArea = None,
                   distanceRange = None
                   ):
        
        (phasename, phaselist) = self._resolve_phasenames(phasesOrPhaseGroupList)
        print("Searching using: %s %s" % (phasename, phaselist), file=sys.stderr)

        if isinstance(networkStationList, str):
            networkStationList = [ networkStationList ]

        if networkStationList is None:
            networkStationList = [ "*.*" ]

        (t0, t1, targetSamplingRate, allowedLocGainList, dataWindowRange, phasesOrPhaseGroupList,
        networkStationList, stationRestrictionArea,
        eventRestrictionArea, magnitudeRange, depthRange,
        distanceRange) = self._check_param(t0, t1, targetSamplingRate, allowedLocGainList, dataWindowRange, phasesOrPhaseGroupList,
                                          networkStationList, stationRestrictionArea,
                                          eventRestrictionArea,
                                          magnitudeRange,
                                          depthRange,
                                          distanceRange)
        
        if self._plotevents:
            print("Cannot plot events yet")
        
        lines = []
        
        for origin, magnitude in self.__loadcsv(t0, t1, eventRestrictionArea, magnitudeRange, depthRange, None, None):
            print("\n Working on origin: %s" % str(origin.time), file=sys.stderr)
            for code in networkStationList:
                (net, sta) = code.split(".")
                
                if self._client_type == "arclink":
                    iterator = self.__loadstArclink(net, sta, t0, t1, stationRestrictionArea, origin, distanceRange)
                elif self._client_type == "fdsn":
                    iterator = self.__loadstFDSN(net, sta, t0, t1, stationRestrictionArea, origin, distanceRange)
                else:
                    raise Exception("Bad Client.")
                
                for network, station in iterator:
                    try:
                        print("  Working on station %s.%s " % (network.code, station.code), end="", file=sys.stderr)
                        self._build(lines,
                                   network, station, origin, magnitude,
                                   phaselist, phasename,
                                   targetSamplingRate, allowedLocGainList, dataWindowRange)
                        print("OK!", file=sys.stderr)
                    except NextItem as e:
                        print("\n  Skipping: %s" % str(e), file=sys.stderr)
        
        request = self._organize_by_event(lines)

        return request
    
    def stationBased(self, t0, t1, targetSamplingRate, allowedLocGainList, dataWindowRange, phasesOrPhaseGroupList,
                    networkStationList,
                    stationRestrictionArea,

                    eventRestrictionArea = None,
                    magnitudeRange = None,
                    depthRange = None,
                    distanceRange = None):

        (t0, t1, targetSamplingRate, allowedLocGainList, dataWindowRange, phasesOrPhaseGroupList,
         networkStationList, stationRestrictionArea,
         eventRestrictionArea, magnitudeRange, depthRange,
         distanceRange) = self._check_param(t0, t1, targetSamplingRate, allowedLocGainList, dataWindowRange, phasesOrPhaseGroupList,
                                          networkStationList, stationRestrictionArea,
                                          eventRestrictionArea,
                                          magnitudeRange,
                                          depthRange,
                                          distanceRange)

        (phasename, phaselist) = self._resolve_phasenames(phasesOrPhaseGroupList)
        print("Searching using: %s %s" % (phasename, phaselist), file=sys.stderr)

        # List of request lines
        lines = []

        # Check if network/station code is a string -> list
        if isinstance(networkStationList, str):
            networkStationList = [ networkStationList ]

        ## Start the Loop
        # On the given station patterns
        for code in networkStationList:
            (net, sta) = code.split(".")

            if self._client_type == "arclink":
                iterator = self.__loadstArclink(net, sta, t0, t1, stationRestrictionArea, None, None)
            else:
                iterator = self.__loadstFDSN(net, sta, t0, t1, stationRestrictionArea, None, None)

            for network, station in iterator:
                print("\n Working on station %s.%s" % (network.code, station.code), file=sys.stderr)
                for origin, magnitude in self.__loadcsv(t0, t1, eventRestrictionArea, magnitudeRange, depthRange, station, distanceRange):
                    try:
                        print("  Working on origin: %s " % str(origin.time), file=sys.stderr, end="")
                        self._build(lines,
                                   network, station, origin, magnitude,
                                   phaselist, phasename,
                                   targetSamplingRate, allowedLocGainList, dataWindowRange)
                        print("OK!", file=sys.stderr)
                    except NextItem as e:
                        print("\n  Skipping: %s" % str(e), file=sys.stderr)

        request = self._organize_by_station(lines)

        return request

'''
  Test Main Code
'''
if __name__ == "__main__":
    from BaseBuilder import Range, AreaRange
#     rb = ArcLinkRequestBuilder("IRIS","seisrequest.iag.usp.br:18001:m.bianchi@iag.usp.br")
#     rb = FDSNBuilder("IRIS")

#     req = rb.load_request("request.pik")

#     Call the stationBased
#     req = rb.eventBased(t0 = UTCDateTime("2007-01-01"),
#                         t1 = UTCDateTime("2008-01-01"),
#                         targetSamplingRate = 20.0,
#                         allowedLocGainList = ["H", "L"],
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
#
#     rb = CSVBuilder("ff", "seismaster:18001:m.bianchi@iag.usp.br")
#     rb = CSVBuilder("ff", "http://seismaster:18003")
#
#     rq = rb.eventBased("2015-01-01", "2018-12-31", 100., ["H"], Range(-200, 200),
#                        "pgroup", AreaRange.WORLD(), Range.ALLMAGS(), Range.ALLDEPTHS(), ["BX.LD*"], None, None)
#
#     rq = rb.stationBased("2015-01-01", "2018-12-31", 100., ["H"], Range(-200, 200), 
#                          "pgroup", [ "BX.LD*" ], AreaRange.WORLD(), None, None, None, None)
#     rb.save_request("teste.req", rq, True)
#
#     rq = rb.load_request("teste.req")
#     BaseBuilder.show_request(rq, True)
    pass
