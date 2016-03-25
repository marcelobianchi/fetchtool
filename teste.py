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

from obspy import UTCDateTime
from RequestBuilder import AreaRange, Range
from RequestBuilder import RequestBuilder, ArcLinkFDSNRequestBuilder
from Downloader import Downloader, FDSNFetcher, Sc3ArclinkFetcher
from Savers import SacSaver, QSaver 
import sys, os

#fetcher = FDSNFetcher("http://www.moho.iag.usp.br/", allinone = False, merge = True)
#fetcher = Sc3ArclinkFetcher("seismaster:18001/oliveira@iag.usp.br", allinone = True, merge = False)
#saver = QSaver(debug = False)
#saver.enableTimeWindowCheck(30.0, 120.0)

#  
# #Call the stationBased
def test_save_load():
    rb = RequestBuilder("IRIS")
    r = rb.eventBased(t0 = UTCDateTime("2014-01-01"),
                      t1 = UTCDateTime("2015-01-01"),
                      targetSamplingRate = 20.0,
                      allowedGainCodes = ["H", "L"],
                      timeRange = Range(-600.0, 600.0),
                      phasesOrPhaseGroup = [ "pgroup" ],
               
                      networkStationCodes = [ "G.ECH"],
                      stationRestrictionArea = AreaRange(-150.0, 60.0, 20.0, 55.0),
               
                      eventRestrictionArea = AreaRange(0.0, 60.0, -30.0, -60.0),
                      magnitudeRange = Range(6.5, 9.0),
                      depthRange = Range(0.0, 400.0),
                      distanceRange = None)
    filename="ev_africa"
    if os.path.isfile(filename):
        os.remove(filename)
    rb.save_request(filename,r)
    r = rb.load_request("%s" % filename)

def test_rb_event():
    rb = RequestBuilder("IRIS")
    r = rb.eventBased(t0 = UTCDateTime("2014-01-01"),
                      t1 = UTCDateTime("2015-01-01"),
                      targetSamplingRate = 20.0,
                      allowedGainCodes = ["H", "L"],
                      timeRange = Range(-600.0, 600.0),
                      phasesOrPhaseGroup = [ "pgroup" ],
               
                      networkStationCodes = [ "G.ECH"],
                      stationRestrictionArea = AreaRange(-150.0, 60.0, 20.0, 55.0),
               
                      eventRestrictionArea = AreaRange(0.0, 60.0, -30.0, -60.0),
                      magnitudeRange = Range(6.5, 9.0),
                      depthRange = Range(0.0, 400.0),
                      distanceRange = None)
    rb.show_request(r)

def test_rb_station():
    rb = RequestBuilder("IRIS")
    r = rb.stationBased(t0 = UTCDateTime("2014-01-01"),
                      t1 = UTCDateTime("2015-01-01"),
                      targetSamplingRate = 20.0,
                      allowedGainCodes = ["H", "L"],
                      timeRange = Range(-600.0, 600.0),
                      phasesOrPhaseGroup = [ "pgroup" ],
              
                      networkStationCodes = [ "G.ECH"],
                      stationRestrictionArea = AreaRange(-150.0, 60.0, 20.0, 55.0),
              
                      eventRestrictionArea = None,
                      magnitudeRange = Range(6.0, 9.0),
                      depthRange = Range(0.0, 400.0),
                      distanceRange = Range(0,90.0))
    rb.show_request(r)

def test_arb_event():
    arb = ArcLinkFDSNRequestBuilder("http://www.moho.iag.usp.br", "seisrequest.iag.usp.br:18001:m.bianchi@iag.usp.br")
    r = arb.eventBased(t0 = UTCDateTime("2014-01-01"),
                      t1 = UTCDateTime("2015-01-01"),
                      targetSamplingRate = 20.0,
                      allowedGainCodes = ["H", "L"],
                      timeRange = Range(-600.0, 600.0),
                      phasesOrPhaseGroup = [ "pgroup" ],
               
                      networkStationCodes = [ "BL.AQDB"],
                      stationRestrictionArea = AreaRange(-150.0, 60.0, 20.0, 55.0),
               
                      eventRestrictionArea = AreaRange(0.0, 60.0, -30.0, -60.0),
                      magnitudeRange = Range(6.0, 9.0),
                      depthRange = Range(0.0, 400.0),
                      distanceRange = None)
    arb.show_request(r)

def test_arb_station():
    arb = ArcLinkFDSNRequestBuilder("http://www.moho.iag.usp.br", "seisrequest.iag.usp.br:18001:m.bianchi@iag.usp.br")
    r = arb.stationBased(t0 = UTCDateTime("2014-01-01"),
                          t1 = UTCDateTime("2015-01-01"),
                          targetSamplingRate = 20.0,
                          allowedGainCodes = ["H", "L"],
                          timeRange = Range(-600.0, 600.0),
                          phasesOrPhaseGroup = [ "pgroup" ],
                  
                          networkStationCodes = [ "BL.AQDB"],
                          stationRestrictionArea = AreaRange(-150.0, 60.0, 20.0, 55.0),
                  
                          eventRestrictionArea = None,
                          magnitudeRange = Range(6.0, 9.0),
                          depthRange = Range(0.0, 400.0),
                          distanceRange = Range(0, 90.0))
    arb.show_request(r)


#===============================================================================
# dl = Downloader(basedir = "/home_ad/thais/IGOR/%s_SacSaver" % filename,
#                  replacetree=False,
#                  resume = True,
#                  fetcher = fetcher,
#                  extracter = [ saver ])
#            
# dl.work(r)
#===============================================================================

if __name__ == "__main__":
    test_rb_event()
    #test_rb_station()
    
    #test_arb_event()
    #test_arb_station()
    
    #test_save_load()
    
    sys.exit()
