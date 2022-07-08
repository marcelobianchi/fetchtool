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

from fetchtool.BaseBuilder import BaseBuilder, NextItem, BadParameter

import sys,socket,os

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
Support Methods - should not be used but can!
'''

def __do_list_networks__(_CL, t0, t1):
        mda = _CL.get_stations(UTCDateTime(t0), UTCDateTime(t1), level = 'net')
        
        for N in mda:
            print("%-02s" % N.code,
                  "[%10s - %10s]" % (N.start_date.strftime("%Y-%m-%d"),N.end_date.strftime("%Y-%m-%d") if N.end_date is not None else "    --    "),
                  " - ", N.description)
        
        return mda


def __do_list_stations__(_CL, t0, t1, net = '*', do_map = False, global_map = False, grids_on = False):
        plot_data = { }
        
        mda = _CL.get_stations(UTCDateTime(t0), UTCDateTime(t1), net = net, level = 'sta')
        
        for N in mda:
            K = "%s_%s" % (N.code, N.start_date.strftime("%Y"))
            plot_data[K] = { 'code' : [], 'lat' : [], 'lon' : []}
            for S in N.stations:
                desc = S.description
                plot_data[K]['lat'].append(S.latitude)
                plot_data[K]['lon'].append(S.longitude)
                plot_data[K]['code'].append(S.code)
                if desc is None:
                    desc = S.site.name if S.site.region is None else "%s (%s)" % (S.site.name, S.site.region)
                print("%02s.%-05s" % (N.code, S.code),
                      "[%10s - %10s]" % (S.start_date.strftime("%Y-%m-%d"),S.end_date.strftime("%Y-%m-%d") if S.end_date is not None else "    --    "),
                      " - ", desc)
            print("")

        if do_map:
            import warnings
            from cartopy import crs, feature as cfeature
            from matplotlib import pyplot as plt
            warnings.filterwarnings('ignore')
            
            fig = plt.figure(figsize=(11, 8.5))
            
            proj = crs.Robinson(central_longitude=0.0)
            data_trans=crs.Geodetic()
            
            ax = plt.subplot(1, 1, 1, projection = proj)
            ax.coastlines()
            ax.add_feature(cfeature.BORDERS, linewidth=0.5, edgecolor='black')
            
            for k in plot_data:
                ax.plot(plot_data[k]['lon'], plot_data[k]['lat'], "^", label = k, transform = data_trans)
            
            if global_map:
                ax.set_global()
            
            if grids_on:
                ax.gridlines(draw_labels=True, dms=True, x_inline=False, y_inline=False)

            ax.legend()
            plt.title('Request Container')
            plt.show()

        return mda


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

    def __str__(self):
        return f'FDSNBuilder w/ events = {self.e_fdsn_client.base_url} and data {self.s_fdsn_client.base_url}'

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

    def __collect_stations(self, E, O):
        
        PA = []
        
        for A in O.arrivals:
            if A.time_weight == 0: continue
            P = [ P for P in E.picks if P.resource_id == A.pick_id ][0]
            PA.append((A, P))
        
        sts = { }
        bulk = []
        
        for A, P in PA:
            N = P.waveform_id.network_code
            S = P.waveform_id.station_code
            L = P.waveform_id.location_code if P.waveform_id.location_code is not None else ""
            C = P.waveform_id.channel_code[:2]
            NS = f'{N}.{S}'
            if NS not in sts:
                sts[NS] = {}
                bulk.append((N,S,L,C + "*", O.time - 3600, O.time + 7200))
            
            if A.phase in sts[NS]:
                raise Exception(f"Duplicate {A.phase} pick for station {NS}")
            
            sts[NS][A.phase] = (N, S, L, C, A.distance, A.azimuth, A.phase, P.time)
        
        return sts, bulk

    '''
    Quick MDA query methods
    '''
    def list_networks(self, t0, t1):
        '''
        Perform a simple query to find networks in the station server
        
        Returns: 
            Metadata object.
        '''
        return __do_list_networks__(self.s_fdsn_client, t0, t1)

    def list_stations(self, t0, t1, net = '*', do_map = False, global_map = False, grids_on = False):
        '''
        Perform a simple query to find stations in the station server. You can
        supply many different networks using wildcards or using the "," as item
        separators.
        
        Returns: 
            Metadata object.
            Optionaly plot the station locations.
        '''
        return __do_list_stations__(self.s_fdsn_client, t0, t1, net, do_map, global_map, grids_on)

    '''
    Request Builder Methods
    '''
    def eventidBased(self, eventid_or_list_of, dataWindowRange):
        lines = []
        
        if isinstance(eventid_or_list_of, str):
            if "," in eventid_or_list_of:
                eventid_or_list_of = eventid_or_list_of.split(",")
            else:
                eventid_or_list_of = [ eventid_or_list_of ]
        
        events = None
        for eid in eventid_or_list_of:
            try:
                kwargsevent = {
                    'eventid' : eid,
                    'includearrivals' : True
                }
                
                event = self.e_fdsn_client.get_events(**kwargsevent)
                
                if events is None:
                    events = event
                    continue
                
                events.extend(event)
            except fdsn.header.FDSNException:
                print(f"No events found for the given parameters {eid}.", file=sys.stderr)

        if len(events) == 0:
            return None
        
        print(f"Total of {len(events)} event found.", file=sys.stderr)
        if self._plotevents:
            events.plot()

        # Event loop
        for event in events:
            try:
                (origin, magnitude) = self._getOrigin(event)
            except NextItem as e:
                print("Skipping Origin: %s" % str(e), file=sys.stderr)
                continue

            print("Working on origin: %s" % str(origin.time), file=sys.stderr)

            try:
                sts,bulk = self.__collect_stations(event, origin)
                inventory = self.s_fdsn_client.get_stations_bulk(bulk, level = 'channel')

                if inventory is None or len(inventory.networks) == 0:
                    print("\n No station could be found.", file=sys.stderr)
                    continue
            except fdsn.header.FDSNException as e:
                print("\n Failed to get information from event to fetch stations", file = sys.stderr)
                continue

            # Event loop
            used = []
            for network in inventory.networks:
                for station in network.stations:
                    try:
                        self._build_predefined(lines,
                                   network, station, origin, magnitude,
                                   sts[f'{network.code}.{station.code}'], 
                                   dataWindowRange)
                        print("  ++ Station %2s.%-5s was added to request." % (network.code, station.code), file=sys.stderr)
                        used.append("%s.%s" % (network.code, station.code))
                    except NextItem as e:
                        print(" -- Skipping: %s" % str(e), file=sys.stderr)
                    except Exception:
                        print(" -- Skipping: %2s.%-5s, general error occurred." % (network.code, station.code), file=sys.stderr)

            for ns in sts:
                if ns not in used:
                    print("  -- Station %-8s was not found the FDSN, skipping." % ns, file=sys.stderr)

        request = self._organize_by_event(lines)

        return request

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
                    print("\nWorking on station %s.%s" % (network.code, station.code), file=sys.stderr)

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
                                # ~ print("Getting events ..")
                                events = self.e_fdsn_client.get_events(**kwargsevent)
                            except socket.timeout:
                                print('Failed to get ev. try %d' % tryid)
                                tryid = tryid + 1
                                pass
                            except socket.error:
                                print('Failed to get ev. try %d' % tryid)
                                tryid = tryid + 1
                                pass
                        if events == "INVALID": events = None
                    except fdsn.header.FDSNException as e:
                        print("FAILED:", e)
                        events = None

                    # Event loop
                    if events is None or len(events) < 0:
                        print(" No Events Found.", file=sys.stderr)
                        continue

                    for event in events:
                        try:
                            (origin, magnitude) = self._getOrigin(event)

                            print(" Working on origin: %s" % str(origin.time), end="", file=sys.stderr)
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

                print(" Stations for pattern: %s.%s" % (net,sta), file=sys.stderr)
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
                            print(" -- Skipping: %s" % str(e), file=sys.stderr)

        request = self._organize_by_event(lines)

        return request


class CSVBuilder(BaseBuilder):
    '''
    Build a request from a Text File for Event and a FDSN/ArcLink server for stations.
    File format is:
    
    <Origin Time> <Longitude> <Latitude> <Depth> <Magnitude> [...]
    
    Separators and also be [ ], [,] or even [;].
    '''

    def __init__(self, time_lon_lat_dep_in_km_mag_file, fdsn_station_url, debug = False):
        BaseBuilder.__init__(self)
        
        self._evfile = time_lon_lat_dep_in_km_mag_file
        self._client = fdsn.Client(fdsn_station_url, debug = debug)
        
        return

    '''
    Quick MDA query methods
    '''
    def list_networks(self, t0, t1):
        '''
        Perform a simple query to find networks in the station server
        
        Returns: 
            Metadata object.
        '''
        return __do_list_networks__(self._client, t0, t1)

    def list_stations(self, t0, t1, net = '*', do_map = False, global_map = False, grids_on = False):
        '''
        Perform a simple query to find stations in the station server. You can
        supply many different networks using wildcards or using the "," as item
        separators.
        
        Returns: 
            Metadata object.
            Optionaly plot the station locations.
        '''
        return __do_list_stations__(self._client, t0, t1, net, do_map, global_map, grids_on)

    '''
    CSV specific Methods
    '''
    def _getChannelList(self, station, t0, targetsps, instcode):
        return FDSNBuilder._getChannelList(station, t0, targetsps, instcode)

    def __loadcsv(self, t0, t1, evr, mr, dr, station = None, distanceRange = None):
        event     = AttribDict()
        magnitude = AttribDict()
        
        if station is not None:
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
                event.depth = dep * 1000.0
                magnitude.mag = mag
                
                '''
                Yield Data
                '''
                yield (event, magnitude)

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
            return
        
        for network in inventory.networks:
            for station in network.stations:
                '''
                Yield Data
                '''
                yield network, station

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
                
                iterator = self.__loadstFDSN(net, sta, t0, t1, stationRestrictionArea, origin, distanceRange)
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
    if not os.path.isfile("catalog.xy"):
        print("Please run prepare_tests.py first!")
        sys.exit(1)

    from BaseBuilder import Range, AreaRange
    
    iris = FDSNBuilder("IRIS")
    print("iris = ", iris)

    usp  = FDSNBuilder("IRIS", "USP")
    print("usp  = ", usp)

    print("\n Using USP server \n")

    print("\n** 3 events with stations from BL.AQDB only\n")
    rq = usp.stationBased("2020-01-01", "2021-01-01", 100, [ "H" ], Range(5,10), "ttp",
                    "BL.AQDB",
                    AreaRange.WORLD(),

                    eventRestrictionArea = AreaRange(0,180,0,90),
                    magnitudeRange = Range(7,10),
                    depthRange = None,
                    distanceRange = None)
    
    usp.show_request(rq)
    
    print("\n** Same 3 events with stations from BL only\n")
    rq = usp.eventBased("2020-01-01", "2021-01-01", 100, [ "H" ], Range(5,10), "ttp",
                   AreaRange(0,180,0,90),
                   Range(7,10),
                   Range.ALLDEPTHS(),

                   networkStationList = [ "BL.*", "BR.*" ],
                   stationRestrictionArea = None,
                   distanceRange = None)
    
    usp.show_request(rq)

    print("\n** Testing CSV class \n")

    csv = CSVBuilder("catalog.xy", "USP")

    rq = csv.stationBased("2020-01-01", "2023-01-01", 100, [ "H" ], Range(5,10), "ttp",
                    "BL.AQDB",
                    AreaRange.WORLD(),

                    eventRestrictionArea = None,
                    magnitudeRange = Range(6,10),
                    depthRange = None,
                    distanceRange = None)

    csv.show_request(rq)

    rq = csv.eventBased("2020-01-01", "2023-01-01", 100, [ "H" ], Range(5,10), "ttp",
                   AreaRange.WORLD(),
                   Range.ALLMAGS(),
                   Range.ALLDEPTHS(),

                   networkStationList = [ "BL.A*" ],
                   stationRestrictionArea = None,
                   distanceRange = Range(0,120))

    csv.show_request(rq)

    rq = csv.stationBased("2020-01-01", "2023-01-01", 100, [ "H" ], Range(5,10), "ttp",
                    "BL.AQDB",
                    AreaRange.WORLD(),

                    eventRestrictionArea = None,
                    magnitudeRange = None,
                    depthRange = Range(0,40),
                    distanceRange = None)

    csv.show_request(rq)
    
