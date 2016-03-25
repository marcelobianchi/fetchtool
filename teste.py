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
from Savers import SacSaver, QSaver , MSSaver
import os

#fetcher = FDSNFetcher("http://www.moho.iag.usp.br/", allinone = False, merge = True)
#fetcher = Sc3ArclinkFetcher("seismaster:18001/oliveira@iag.usp.br", allinone = True, merge = False)
#saver = QSaver(debug = False)
#saver.enableTimeWindowCheck(30.0, 120.0)

def test_save_load():
    print(" ** Test Request Builder **")
    print(" ** Test Load - Save     **")
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
    print ("\n\n")
    return r

def test_rb_event():
    print(" ** Test Request Builder **")
    print(" ** Event Based **")
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
    print ("\n\n")
    return r

def test_rb_station():
    print(" ** Test Request Builder **")
    print(" ** Station Based **")
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
                      magnitudeRange = Range(7.0, 9.0),
                      depthRange = Range(0.0, 400.0),
                      distanceRange = Range(0,90.0))
    print ("\n\n")
    return r

def test_arb_event():
    print(" ** Test ArcLink Request Builder **")
    print(" ** Event Based **")
    arb = ArcLinkFDSNRequestBuilder("http://www.moho.iag.usp.br", "seisrequest.iag.usp.br:18001:m.bianchi@iag.usp.br")
    r = arb.eventBased(t0 = UTCDateTime("2014-01-01"),
                      t1 = UTCDateTime("2015-01-01"),
                      targetSamplingRate = 20.0,
                      allowedGainCodes = ["H", "L"],
                      timeRange = Range(-600.0, 600.0),
                      phasesOrPhaseGroup = [ "pgroup" ],
               
                      networkStationCodes = [ "BL.AQDB", "BL.PCMB"],
                      stationRestrictionArea = AreaRange(-150.0, 60.0, 20.0, 55.0),
               
                      eventRestrictionArea = AreaRange(0.0, 60.0, -30.0, -60.0),
                      magnitudeRange = Range(6.0, 9.0),
                      depthRange = Range(0.0, 400.0),
                      distanceRange = None)
    print ("\n\n")
    return r

def test_arb_station():
    print(" ** Test ArcLink Request Builder **")
    print(" ** Station Based **")
    arb = ArcLinkFDSNRequestBuilder("http://www.moho.iag.usp.br", "seisrequest.iag.usp.br:18001:m.bianchi@iag.usp.br")
    r = arb.stationBased(t0 = UTCDateTime("2014-01-01"),
                          t1 = UTCDateTime("2015-01-01"),
                          targetSamplingRate = 20.0,
                          allowedGainCodes = ["H", "L"],
                          timeRange = Range(-600.0, 600.0),
                          phasesOrPhaseGroup = [ "pgroup" ],
                  
                          networkStationCodes = [ "BL.AQDB", "BL.PCMB" ],
                          stationRestrictionArea = AreaRange(-150.0, 60.0, 20.0, 55.0),
                  
                          eventRestrictionArea = None,
                          magnitudeRange = Range(7.0, 9.0),
                          depthRange = Range(0.0, 400.0),
                          distanceRange = Range(0, 90.0))
    print ("\n\n")
    return r

if __name__ == "__main__":
    test_save_load()
    
    r = test_rb_station()
    r = test_rb_event()
    
    rs = test_arb_station()
    re = test_arb_event()
    
    print(" ** Test Request Builder **")
    print(" ** Show Request **")
    ArcLinkFDSNRequestBuilder.show_request(r, False)
    
    print(" ** Test Request Builder **")
    print(" ** Show Request Compact **")
    ArcLinkFDSNRequestBuilder.show_request(r)
    
    print(" ** Test Fetching **")
    print(" ** FDSN fetcher **")
    fetcher = FDSNFetcher("http://seisrequest.iag.usp.br/",
                          allinone = False,
                          merge = True)
    
    dl = Downloader(basedir = "./Test_Save/",
                      replacetree=True,
                      show_resume = False,
                      fetcher = fetcher,
                      extracter = None)
    dl.work(r)
    print("\n\n")

    print(" ** Test Fetching **")
    print(" ** ArcLink fetcher **")
    fetcher = Sc3ArclinkFetcher("seisrequest.iag.usp.br:18001/m.bianchi@iag.usp.br",
                          allinone = False,
                          merge = True)
    
    dl = Downloader(basedir = "./Test_Save/",
                      replacetree=True,
                      show_resume = False,
                      fetcher = fetcher,
                      extracter = None)
    dl.work(r)
    print("\n\n")
    
    print(" ** Test Saving (Sac + Qfile + MS modes 1,2,3) **")
    print(" ** ArcLink fetcher **")
    
    assac = SacSaver()
    asq   = QSaver()
    
    asm1   = MSSaver(1) # Mode 1
    asm2   = MSSaver(2) # Mode 2
    asm3   = MSSaver(3) # Mode 3
    
    dl = Downloader(basedir = "./Test_SaveE/",
                      replacetree=True,
                      show_resume = True,
                      fetcher = fetcher,
                      extracter = [ assac, asq, asm1, asm2, asm3 ])
    dl.enableSaveRaw()
    dl.work(re)
    
    dl = Downloader(basedir = "./Test_SaveS/",
                      replacetree=True,
                      show_resume = True,
                      fetcher = fetcher,
                      extracter = [ assac, asq, asm1, asm2, asm3 ])
    dl.work(rs)
