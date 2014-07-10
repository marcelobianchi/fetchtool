from obspy import UTCDateTime
from RequestBuilder import AreaRange, Range
from RequestBuilder import RequestBuilder
from Downloader import Downloader, FDSNFetcher
from Savers import SacSaver, QSaver 
import sys

# Get an Instance of the class
rb = RequestBuilder("IRIS")
fetcher = FDSNFetcher("IRIS", allinone = False, merge = True)
saver = QSaver(debug = False)
saver.enableTimeWindowCheck(350.0, 200.0)

filename="request-adquirindo_n_dados-6.5-9.0_2007_2014"
# 
# #Call the stationBased
# r = rb.eventBased(t0 = UTCDateTime("2007-01-01"),
#                  t1 = UTCDateTime("2015-01-01"),
#                  targetSamplingRate = 20.0,
#                  allowedGainCodes = ["H", "L"],
#                  timeRange = Range(-600.0, 600.0),
#                  phasesOrPhaseGroup = [ "SS" ],
#       
#                  networkStationCodes = [ "TA.*"],
#                  stationRestrictionArea = AreaRange(-150.0, -90.0, 15.0, 60.0),
#       
#                  eventRestrictionArea = AreaRange(-35.0, -10.0, -55.0, -60.0),
#                  magnitudeRange = Range(6.5, 9.0),
#                  depthRange = Range(0.0, 400.0),
#                  distanceRange = None,
#       
#                  )
# rb.save_request(filename,r)
# # sys.exit()

r = rb.load_request("%s" % filename)

for key in r:
   for l in r[key]:
       print l[4]
# print r.keys()
# for i in ['140629_075256', '100105_045538', '110306_143236', '080630_061743', '080223_155719']:
#      del r[i]
# print r.keys()
#   
# dl = Downloader(basedir = "/home_ad/thais/SS/%s_Qsaver" % filename,
#                 replacetree=True,
#                 resume = True,
#                 fetcher = fetcher,
#                 extracter = [ saver ])
#        
# dl.work(r)
sys.exit()
