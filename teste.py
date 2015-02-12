from obspy import UTCDateTime
from RequestBuilder import AreaRange, Range
from RequestBuilder import RequestBuilder
from Downloader import Downloader, FDSNFetcher, Sc3ArclinkFetcher
from Savers import SacSaver, QSaver 
import sys

# Get an Instance of the class
rb = RequestBuilder("IRIS")
fetcher = FDSNFetcher("http://www.moho.iag.usp.br/", allinone = False, merge = True)
#fetcher = Sc3ArclinkFetcher("seismaster:18001/oliveira@iag.usp.br", allinone = True, merge = False)
saver = QSaver(debug = False)

saver.enableTimeWindowCheck(30.0, 120.0)

filename="ev_africa"
#  
# #Call the stationBased
# r = rb.eventBased(t0 = UTCDateTime("2007-01-01"),
#                   t1 = UTCDateTime("2015-01-01"),
#                   targetSamplingRate = 20.0,
#                   allowedGainCodes = ["H", "L"],
#                   timeRange = Range(-600.0, 600.0),
#                   phasesOrPhaseGroup = [ "SS" ],
#           
#                   networkStationCodes = [ "TA.*"],
#                   stationRestrictionArea = AreaRange(-150.0, 60.0, 20.0, 55.0),
#           
#                   eventRestrictionArea = AreaRange(0.0, 60.0, -30.0, -60.0),
#                   magnitudeRange = Range(6.5, 9.0),
#                   depthRange = Range(0.0, 400.0),
#                   distanceRange = None,
#           
#                   )
# rb.save_request(filename,r)
# # sys.exit()
 
r = rb.load_request("%s" % filename)
   
for key in r:
    for l in r[key]:
        print l[6]
print r.keys()
# for i in ['110306_143236', '140629_075256']:
#      del r[i]
# print r.keys()
    
#===============================================================================
# dl = Downloader(basedir = "/home_ad/thais/IGOR/%s_SacSaver" % filename,
#                  replacetree=False,
#                  resume = True,
#                  fetcher = fetcher,
#                  extracter = [ saver ])
#            
# dl.work(r)
#===============================================================================
sys.exit()
