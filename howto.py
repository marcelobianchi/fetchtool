from RequestBuilder import RequestBuilder
from RequestBuilder import ArcLinkRequestBuilder
from RequestBuilder import Range, AreaRange
from obspy.core.utcdatetime import UTCDateTime

from Savers import SacSaver
from Downloader import Downloader, Sc3ArclinkFetcher
import sys

rb = ArcLinkRequestBuilder("http://www.moho.iag.usp.br", "rsis1.on.br:18001:m.bianchi@iag.usp.br")

# request = rb.eventBased(t0 = UTCDateTime("2014-03-16"),
#                         t1 = UTCDateTime("2014-05-01"),
#                         targetSamplingRate = 50.0,
#                         allowedGainCodes = [ "H", "L" ],
#                         timeRange = Range(-120.0, 35*60),
#                         phasesOrPhaseGroup = "pgroup",
# 
#                         eventRestrictionArea = AreaRange(-80, -70, -22, -19),
#                         magnitudeRange = None,
#                         depthRange = None,
# 
#                         networkStationCodes = [ "ON.*" ],
#                         stationRestrictionArea = None,
#                         distanceRange = None
#                         )
# rb.save_request("req-on.pik", request)
# sys.exit()

request = rb.load_request("req-on.pik")
rb.show_request(request)

ft = Sc3ArclinkFetcher("rsis1.on.br:18001/m.bianchi@iag.usp.br", allinone = True, merge = False)
sv = SacSaver(debug = False)

params = sv.getDefaultParameters()
params.tw.prephasevalue = 60.0
params.tw.postphasevalue = 600.0

dl = Downloader("/home/mbianchi/ON-CHILE", replacetree = True, resume = True, fetcher = ft, extracter = sv)
dl.work(request, params)
