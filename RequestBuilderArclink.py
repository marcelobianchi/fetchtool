## Python Basics
import sys, pickle

from RequestBuilder import BaseBuilder, Range, AreaRange, NextItem

## Obspy tools
from obspy.core import UTCDateTime, util, AttribDict, read as oREAD, Stream
from obspy import fdsn
from obspy.fdsn.header import FDSNException

## Sc3 Tools
from seiscomp.arclink.manager import ArclinkManager, ArclinkError
from seiscomp.arclink.client import arclink_status_string

# Defaults Time Constants
SECOND = 1
MINUTE = 60*SECOND
HOUR = 60*MINUTE
DAY = 24*HOUR
WEEK = 7 * DAY
MONTH = 30 * DAY
HEAR = 365 * DAY

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

class ArcLinkRequestBuilder(BaseBuilder):
    def __init__(self, fdsnURL, arclinkURL):
        (host,port,user) = arclinkURL.strip().split(":")
        self.arclink_manager = ArclinkManager("%s:%s" % (host,port), user)
        self.fdsn_client = fdsn.Client(fdsnURL)

    def __chooseArcLink(self, item, location, channel, targetsps, instcode):

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
        elif newitem[5] < item[5] and newitem[5] > targetsps:
            return newitem

        return item

    def __getArcLinkChannelList(self, station, t0, targetsps, instcode):
        clist = [ ]

        ## Choose will build something like: (loca, chan, sps , az  , dip , dsps)
        z = (None, None, None, None, None, None)
        n = (None, None, None, None, None, None)
        e = (None, None, None, None, None, None)

        for (lcode,lstart,location) in unWrapNSLC(station.sensorLocation):
            for (ccode,cstart,channel) in unWrapNSLC(location.stream):
                if t0 < channel.start or (channel.end is not None and t0 > channel.end):
                    continue
    
                if channel.code[-1] == "Z":
                    z = self.__chooseArcLink(z, location, channel, targetsps, instcode)
                elif channel.code[-1] == "N" or channel.code[-1] == "1":
                    n = self.__chooseArcLink(n, location, channel, targetsps, instcode)
                elif channel.code[-1] == "E" or channel.code[-1] == "2":
                    e = self.__chooseArcLink(e, location, channel, targetsps, instcode)
    #             else:
    #                 print >>sys.stderr, "Unknow channel %s.%s" % (channel.location_code, channel.code)

        if z[1] is None or n[1] is None or e[1] is None:
            raise NextItem("No Z, N or E channel(s) found")

        # Return only (location, channel, azimuth, dip)
        z = (z[0], z[1], z[3], z[4])
        n = (n[0], n[1], n[3], n[4])
        e = (e[0], e[1], e[3], e[4])

        return (z,n,e)

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
            print >>sys.stderr,""
            origin = self.getOrigin(event)

            print >>sys.stderr,"Working on origin: %s" % str(origin.time)

            for code in networkStationCodes:
                (net, sta) = code.split(".")

                inventory = self.arclink_manager.get_inventory(net, sta, "*", "*", (origin.time - 86400).datetime, (origin.time + 86400).datetime)
                if inventory is None or len(inventory.network) == 0:
                    print >>sys.stderr,"No stations for pattern: %s.%s" % (net, sta)
                    continue

                print >>sys.stderr,"\n Stations for pattern: %s.%s" % (net,sta)
                # Station loop
                for (ncode, nstart, network) in  unWrapNSLC(inventory.network):
                    for (scode, sstart, station) in  unWrapNSLC(network.station):
                        try:
                            print >>sys.stderr,"  Working on station %s.%s " % (network.code, station.code),

                            delta = util.locations2degrees(station.latitude,
                                                           station.longitude,
                                                           origin.latitude,
                                                           origin.longitude)

                            (ta, slowness) = self.find_arrivalTime(origin.time, delta, origin.depth, phaselist)

                            (z,n,e) = self.__getArcLinkChannelList(station, ta, targetSamplingRate, allowedGainCodes)

                            EI = self.build_event_dictionary(origin.time, origin.latitude, origin.longitude, origin.depth)
                            SI = self.build_station_dictionary(network.code, station.code, station.latitude, station.longitude, station.elevation)
                            PI = self.build_pick_dictionary(phasename, ta, slowness)

                            lines.append((ta + timeRange.min(),
                                          ta + timeRange.max(),
                                          network.code,
                                          station.code,
                                          [z,n,e],
                                          SI,
                                          EI,
                                          PI
                                         )
                            )
                            print >>sys.stderr,"OK!"
                        except NextItem, e:
                            print >>sys.stderr,"\n  Skipping: %s" % str(e)

        request = self.organize_by_event(lines)

        return request

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

        ## Start the Loop
        # On the given station patterns
        for code in networkStationCodes:
            (net, sta) = code.split(".")
            inventory = self.arclink_manager.get_inventory(net, sta, "*", "*", t0.datetime, t1.datetime)

            if inventory is None or len(inventory.network) == 0:
                print >>sys.stderr,"No stations for pattern: %s.%s" % (net, sta)
                continue

            print >>sys.stderr,"Stations for pattern: %s.%s" % (net,sta)
            # Station Loop
            for (ncode, nstart, network) in  unWrapNSLC(inventory.network):
                for (scode, sstart, station) in  unWrapNSLC(network.station):
                    print >>sys.stderr,"\n Working on station %s.%s" % (network.code, station.code)

                    # Start Build the parameters to pass on the Event Client
                    start = UTCDateTime(station.start) if station.start else t0
                    end = UTCDateTime(station.end) if station.end else t1
                    kwargsevent = self.fill_kwargsevent(max([start, t0]),
                                                        min([end, t1]),
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

                    if events is None or len(events) < 0:
                        print >>sys.stderr,"  No Events Found."
                        continue

                    # Event loop
                    for event in events:
                        try:
                            origin = self.getOrigin(event)

                            print >>sys.stderr,"  Working on origin: %s" % str(origin.time),

                            delta = util.locations2degrees(station.latitude,
                                                           station.longitude,
                                                           origin.latitude,
                                                           origin.longitude)

                            (ta, slowness) = self.find_arrivalTime(origin.time, delta, origin.depth, phaselist)

                            (z,n,e) = self.__getArcLinkChannelList(station, ta, targetSamplingRate, allowedGainCodes)

                            EI = self.build_event_dictionary(origin.time, origin.latitude, origin.longitude, origin.depth)
                            SI = self.build_station_dictionary(network.code, station.code, station.latitude, station.longitude, station.elevation)
                            PI = self.build_pick_dictionary(phasename, ta, slowness)

                            lines.append((ta + timeRange.min(),
                                      ta + timeRange.max(),
                                      network.code,
                                      station.code,
                                      [z,n,e],
                                      SI,
                                      EI,
                                      PI)
                                     )
                            print >>sys.stderr,"OK!"
                        except NextItem,e:
                            print >>sys.stderr,"  Skipping: %s" % str(e)
                            continue

        request = self.organize_by_station(lines)

        return request

if __name__ == "__main__":
    # Get an Instance of the class
    rb = ArcLinkRequestBuilder("IRIS","seisrequest.iag.usp.br:18001:m.bianchi@iag.usp.br")

    # Call the stationBased
    req = rb.eventBased(t0 = UTCDateTime("2007-01-01"),
                    t1 = UTCDateTime("2008-01-01"),
                    targetSamplingRate = 20.0,
                    allowedGainCodes = ["H", "L"],
                    timeRange = Range(-120, 600),
                    phasesOrPhaseGroup = "pgroup",

                    networkStationCodes = [ "BL.*", "TA.*"],
                    stationRestrictionArea = AreaRange(-150.0, -90.0, 15.0, 60.0),

                    eventRestrictionArea = AreaRange(-35.0, -10.0, -55.0, -60.0),
                    magnitudeRange = Range(6., 9.0),
                    depthRange = Range(0.0, 400.0),
                    distanceRange = None
                    )
#     rb.show_request(req)

