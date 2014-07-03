## Python Basics
import sys, pickle

from RequestBuilder import BaseBuilder, Range, AreaRange, NextItem

## Obspy tools
from obspy.core import UTCDateTime, util, AttribDict
from obspy import fdsn
from obspy.station import Station as FDSNStation
from obspy.fdsn.header import FDSNException

## Sc3 Tools
from seiscomp.arclink.manager import ArclinkManager
import obspy

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
        elif newitem[5] < item[5] and newitem[5] > targetsps:
            return newitem

        return item

    def getChannelList(self, station, t0, targetsps, instcode):
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
                    z = self.__choose(z, location, channel, targetsps, instcode)
                elif channel.code[-1] == "N" or channel.code[-1] == "1":
                    n = self.__choose(n, location, channel, targetsps, instcode)
                elif channel.code[-1] == "E" or channel.code[-1] == "2":
                    e = self.__choose(e, location, channel, targetsps, instcode)
    #             else:
    #                 print >>sys.stderr, "Unknow channel %s.%s" % (channel.location_code, channel.code)

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

    def eventBased(self, t0, t1, targetSamplingRate, allowedGainCodes, timeRange, phasesOrPhaseGroup,
                   eventRestrictionArea,
                   magnitudeRange,
                   depthRange,

                   networkStationCodes = None,
                   stationRestrictionArea = None,
                   distanceRange = None
                   ):

        # List of request lines
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
                            self.build(lines,
                                         network, station, origin,
                                         phaselist, phasename,
                                         targetSamplingRate, allowedGainCodes, timeRange)
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

        # List of request lines
        lines = []

        (phasename, phaselist) = self.resolve_phasenames(phasesOrPhaseGroup)
        print >>sys.stderr,"Searching using: %s %s" % (phasename, phaselist)

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

if __name__ == "__main__":
    # Get an Instance of the class
    rb = ArcLinkRequestBuilder("IRIS","seisrequest.iag.usp.br:18001:m.bianchi@iag.usp.br")

    # Call the stationBased
    req = rb.stationBased(t0 = UTCDateTime("2007-01-01"),
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

