from __future__ import division, print_function
from fetchtool.BaseBuilder import BaseBuilder, NextItem, BadParameter
import sys,socket

from seiscomp.arclink.manager import ArclinkManager
from obspy.clients import fdsn

'''
  Tools
'''
from obspy import UTCDateTime
from obspy.core import AttribDict
from obspy import geodetics


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


'''
  Builders Definition
'''
class ArcLinkFDSNBuilder(BaseBuilder):
    '''
    Build a request using a FDSN server for events and an ArcLink for stations metadata
    '''
    def __init__(self, fdsn_event_url, arclink_url):
        BaseBuilder.__init__(self)
        (host,port,user) = arclink_url.strip().replace("/",":").split(":")
        self.__arclink_manager = ArclinkManager("%s:%s" % (host,port), user)
        self.__fdsn_client = fdsn.Client(fdsn_event_url)

    '''
    ArcLink specific Methods
    '''
    @staticmethod
    def _choose(item, location, channel, targetsps, instcode):
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
        sps = channel.sampleRateNumerator / channel.sampleRateDenominator
        newitem = (location.code, channel.code, sps, channel.azimuth, channel.dip, abs(sps - targetsps))

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

        for (_, _, location) in unWrapNSLC(station.sensorLocation):
            for (_, _, channel) in unWrapNSLC(location.stream):
                if t0 < channel.start or (channel.end is not None and t0 > channel.end):
                    continue

                if channel.code[-1] == "Z":
                    z = ArcLinkFDSNBuilder._choose(z, location, channel, targetsps, instcode)
                elif channel.code[-1] == "N" or channel.code[-1] == "1":
                    n = ArcLinkFDSNBuilder._choose(n, location, channel, targetsps, instcode)
                elif channel.code[-1] == "E" or channel.code[-1] == "2":
                    e = ArcLinkFDSNBuilder._choose(e, location, channel, targetsps, instcode)
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

    def eventBased(self, t0, t1, targetSamplingRate, allowedLocGainList, dataWindowRange, phasesOrPhaseGroupList,
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
            except NextItem:
                print("Skipping Origin", file=sys.stderr)
                continue

            print("\nWorking on origin: %s" % str(origin.time), file=sys.stderr)

            for code in networkStationList:
                (net, sta) = code.split(".")
                inventory = self.__arclink_manager.get_inventory(net, sta, "*", "*", (origin.time - 86400).datetime, (origin.time + 86400).datetime)

                if inventory is None or len(inventory.network) == 0:
                    print("No stations for pattern: %s.%s" % (net, sta), file=sys.stderr)
                    continue

                print(" Stations for pattern: %s.%s" % (net,sta), file=sys.stderr)
                # Station loop
                for (_, _, network) in  unWrapNSLC(inventory.network):
                    for (_, _, station) in  unWrapNSLC(network.station):
                        # ! MISSING CHECKS
                        try:
                            print("  Working on station %s.%s " % (network.code, station.code), end="", file=sys.stderr)
                            self._build(lines,
                                         network, station, origin, magnitude,
                                         phaselist, phasename,
                                         targetSamplingRate, allowedLocGainList, dataWindowRange)
                            print("OK!", file=sys.stderr)
                        except NextItem:
                            print("  Skipping.", file=sys.stderr)

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

        # List of request lines
        lines = []

        (phasename, phaselist) = self._resolve_phasenames(phasesOrPhaseGroupList)
        print("Searching using: %s %s" % (phasename, phaselist), file=sys.stderr)

        (t0, t1, targetSamplingRate, allowedLocGainList, dataWindowRange, phasesOrPhaseGroupList,
        networkStationList, stationRestrictionArea,
        eventRestrictionArea, magnitudeRange, depthRange,
        distanceRange) = self._check_param(t0, t1, targetSamplingRate, allowedLocGainList, dataWindowRange, phasesOrPhaseGroupList,
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

            print("\nStations for pattern: %s.%s" % (net,sta), file=sys.stderr)
            # Station Loop
            for (_, _, network) in  unWrapNSLC(inventory.network):
                for (_, _, station) in  unWrapNSLC(network.station):
                    print(" Working on station %s.%s" % (network.code, station.code), file=sys.stderr)

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
                    except fdsn.header.FDSNException:
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
                                       targetSamplingRate, allowedLocGainList, dataWindowRange)
                            print("OK!", file=sys.stderr)
                        except NextItem:
                            print("  Skipping.", file=sys.stderr)
                            continue

        request = self._organize_by_station(lines)

        return request

