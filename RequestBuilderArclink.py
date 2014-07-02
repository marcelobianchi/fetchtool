## Python Basics
import sys, pickle

from RequestBuilder import BaseBuilder, Range, AreaRange

## Obspy tools
from obspy.core import UTCDateTime, util, AttribDict, read as oREAD, Stream
from obspy import fdsn, taup, arclink, sh

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
        self._am = ArclinkManager("%s:%s" % (host,port), user)
        self._fc = fdsn.Client(fdsnURL)

    '''
    To clean methods
    '''
    def _loadInventory(self):
        print >>sys.stderr,"Loading Complete Inventory ..."
        inv = self.am.get_inventory("*", "*", "*", "*", t0.datetime, t1.datetime)
        print >>sys.stderr," filtering inventory ..."
        self._filter(inv)
        print >>sys.stderr,"Inventory if ready.\n"
        return inv

    def _preparePattern(self, pattern):
        return "^" + pattern.replace(".","[.]").replace("*",".*").replace("?",".") + "$"

    def _filter(self, inv):
        if self.selectors == None: return

        for (ncode, nstart, net) in  unWrapNSLC(inv.network):
            for (scode, sstart, sta) in  unWrapNSLC(net.station):
                code = "%s.%s" % (ncode,scode)
                for pattern in self.selectors:
                    if re.match(self._preparePattern(pattern), code) is not None:
                        break
                else:
                    net.remove_station(scode, sstart)

    def __ArclinkChoose(self, item, location, channel, targetsps, instcode):

        newccode = channel.code
        newsps = (channel.sampleRateNumerator / channel.sampleRateDenominator)
        newitem = ( (channel.code, channel.start, channel.azimuth, channel.dip), (location.code, location.start), abs(targetsps - newsps))

        # Check that the instrument code is allowed to select
        if channel.code[1] not in instcode:
            return item

        # If we did not select yet return the current
        if item[0] == None:
            return newitem

        index = item[2]
        ccode  = item[0][0]

        # Decided based on SPS first
        if abs(targetsps - newsps) == index:
            # And based on the Channel Instrument Code
            if instcode.index(newccode[1]) < instcode.index(ccode[1]):
                return newitem
        elif abs(targetsps - newsps) < index:
            return newitem

        return item

    def __getArclinkChannelList(self, station, tp, targetsps, instcode):
        clist = [ ]

        z = (None, None, None)
        n = (None, None, None)
        e = (None, None, None)

        for (lcode, lstart, location) in unWrapNSLC(station.sensorLocation):
            for (ccode, cstart, channel) in unWrapNSLC(location.stream):

                # Check time range
                if (tp.datetime < channel.start) or (channel.end is not None and tp.datetime > channel.end):
                    continue

                # Enter selection of channel
                if ccode[-1] == "Z":
                    z = self._choose(z, location, channel, targetsps, instcode)
                elif ccode[-1] == "N" or ccode[-1] == "1":
                    n = self._choose(n, location, channel, targetsps, instcode)
                elif ccode[-1] == "E" or ccode[-1] == "2":
                    e = self._choose(e, location, channel, targetsps, instcode)
                else:
                    print >>sys.stderr, "Unknow channel %s.%s" % (lcode, ccode)

        if z[0] is not None:
            clist.append((z[1][0], z[0][0], z[0][2], z[0][3]))

        if n[0] is not None:
            clist.append((n[1][0], n[0][0], n[0][2], n[0][3]))

        if e[0] is not None:
            clist.append((e[1][0], e[0][0], e[0][2], e[0][3]))

        return clist

    def _checkEvent(self, event):

        if len(event.origins) == 0:
            print >>sys.stderr,"  Event with no origin"
            return False

        origin = event.preferred_origin()
        if origin is None:
            origin = event.origins[0]

        if origin.depth is None: 
            print >>sys.stderr,"  Origin has no depth"
            return False

        return True

    def getDefaultsParameters(self):
        parameters = AttribDict({})

        # Channel Selections
        parameters['targetsps'] = 20.0
        parameters['prewindow'] = 2.0 * MINUTE
        parameters['postwindow'] = 15.0 * MINUTE
        parameters['instcode'] = ["H", "L"]

        # Event Selection
        parameters['minradius'] = 0.0
        parameters['maxradius'] = 90.0
        parameters['minmagnitude'] = 5.5
        return parameters

    def OLDstationBased(self, parameters):
        request = {}

        print >>sys.stderr,"StationBuilder:"
        for (ncode, nstart, network) in unWrapNSLC(self.inventory.network):
            for (scode, sstart, station) in unWrapNSLC(network.station):
                lines = []
                print >>sys.stderr," Working on %s.%s (" % (ncode,scode),

                # Build Station Info
                stationinfo = self._buildStationInfo(network, station)

                # Fetch Events
                catalog = self._events(station, parameters)
                if not catalog:
                    print >>sys.stderr,"No Events)"
                    continue
                else:
                    print >>sys.stderr,"%d events)" % (len(catalog))

                # Loop on Event
                for event in catalog:
                    if not self._checkEvent(event): 
                        print >>sys.stderr,"  skipping event %s" % event
                        print event.origins
                        continue

                    eventinfo = self._buildEventInfo(event)

                    # Compute arrival time
                    (tp, slowness) = self._arrivalTime(station, event)
                    if tp == None:
                        print >>sys.stderr,"No matching phase for %s on %s" % (event, station.code)

                    # Find matching channel close to SPS
                    channels = self._getChannelList(station, tp, parameters.targetsps, parameters.instcode)
                    if len(channels) == 0:
                        print >>sys.stderr,"  no channels operating at event time = %s" % tp
                        continue

                    # Build pick info
                    pickinfo = self._buildPickInfo("P", tp, slowness)

                    lines.append((tp - parameters.prewindow, tp + parameters.postwindow, ncode, scode, channels, stationinfo, eventinfo, pickinfo))

                if lines:
                    request["%s.%s" % (ncode, scode)] = lines
                else:
                    print >>sys.stderr,"Will not request any data for %s.%s" % (ncode,scode)

        print >>sys.stderr,"StationBuilder is ready.\n"
        ## { STKEY: (T0,T1,N,S,[(LO,C)],At-Station, At-Event, At-Pick) }
        return request

    '''
    New Methods
    '''
    def eventBased(self, t0, t1, targetSamplingRate, allowedGainCodes, timeRange,
                   eventRestrictionArea,
                   magnitudeRange,
                   depthRange,

                   networkStationCodes = None,
                   stationRestrictionArea = None,
                   distanceRange = None
                   ):

        # List of request lines
        lines = []

        kwargsevent = self.fill_kwargsevent(t0, t1,
                                              eventRestrictionArea,
                                              magnitudeRange,
                                              depthRange,
                                              None, None, None)

        try:
            events = self._fc.get_events(**kwargsevent)
            events.plot()
            print >>sys.stderr,"Found %d events." % len(events)
        except fdsn.header.FDSNException:
            print >>sys.stderr,"\nNo events found for the given parameters.\n"
            sys.exit()

        for event in events:
            origin = event.preferred_origin()
            if origin is None:
                print >>sys.stderr,"Bad origin."
                continue

        return

    def stationBased(self, t0, t1, targetSamplingRate, allowedGainCodes, timeRange,
                    networkStationCodes,
                    stationRestrictionArea,

                    eventRestrictionArea = None,
                    magnitudeRange = None,
                    depthRange = None,
                    distanceRange = None):

        # List of request lines
        lines = []
        return

if __name__ == "__main__":
    # Get an Instance of the class
    rb = ArcLinkRequestBuilder("IRIS","seisrequest.iag.usp.br:18001:m.bianchi@iag.usp.br")

    # Call the stationBased
    req = rb.eventBased(t0 = UTCDateTime("2007-01-01"),
                    t1 = UTCDateTime("2008-01-01"),
                    targetSamplingRate = 20.0,
                    allowedGainCodes = ["H", "L"],
                    timeRange = Range(-120, 600),

                    networkStationCodes = [ "TA.*"],
                    stationRestrictionArea = AreaRange(-150.0, -90.0, 15.0, 60.0),

                    eventRestrictionArea = AreaRange(-35.0, -10.0, -55.0, -60.0),
                    magnitudeRange = Range(6.2, 9.0),
                    depthRange = Range(0.0, 400.0),
                    distanceRange = None
                    )


