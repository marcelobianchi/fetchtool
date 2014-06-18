from obspy import UTCDateTime
from RequestBuilder import AreaRange, Range
from RequestBuilder import RequestBuilder
from Downloader import Downloader, FDSNFetcher
from Savers import SacSaver 
import sys

# Get an Instance of the class
rb = RequestBuilder("IRIS")
fetcher = FDSNFetcher("IRIS", allinone = False, merge = True)
saver = SacSaver(debug = False)
saverparam = saver.getDefaultParameters()

saverparam.tw.prephasevalue  = 10.0
saverparam.tw.postphasevalue = 40.0
saverparam.rms = None 
saverparam.rmsratio = 2.0


# Call the stationBased
r = rb.eventBased(t0 = UTCDateTime("2011-01-01"),
                t1 = UTCDateTime("2011-02-01"),
                targetSamplingRate = 20.0,
                allowedGainCodes = ["H", "L"],
                timeRange = Range(-20.0, 50.0),

                networkStationCodes = [ "TA.Z*"],
                stationRestrictionArea = AreaRange(-150.0, -90.0, 15.0, 35.0),

                eventRestrictionArea = AreaRange(-75.0, -15.0, -35.0, -45.0),
                magnitudeRange = Range(5.0, 7.0),
                depthRange = Range(0.0, 400.0),
                distanceRange = None
                )

dl = Downloader(basedir = "/home_ad/mbianchi/datateste/",
                replacetree=True,
                resume = True,
                fetcher = fetcher,
                extracter = [ saver ])

dl.work(r, saverparam)
