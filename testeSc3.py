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
from BaseBuilder import AreaRange, Range
from BuildersSc3 import ArcLinkFDSNBuilder
from Downloader import Downloader
from DownloaderSc3 import Sc3ArclinkFetcher
from Savers import SacSaver, QSaver , MSSaver

import os


if __name__ == "__main__":
    print("\n")
    print(" ** Test ArcLink Request Builder **")
    print(" ** Event Based **")
    print("\n")
    arb = ArcLinkFDSNBuilder("http://www.moho.iag.usp.br", "seisrequest.iag.usp.br:18001:m.bianchi@iag.usp.br")
    rs = arb.eventBased(t0 = UTCDateTime("2014-01-01"),
                       t1 = UTCDateTime("2015-01-01"),
                       targetSamplingRate = 20.0,
                       allowedLocGainList = ["H", "L"],
                       dataWindowRange = Range(-600.0, 600.0),
                       phasesOrPhaseGroupList = [ "pgroup" ],
               
                       networkStationList = [ "BL.AQDB", "BL.PCMB"],
                       stationRestrictionArea = AreaRange(-150.0, 60.0, 20.0, 55.0),
               
                       eventRestrictionArea = AreaRange(0.0, 60.0, -30.0, -60.0),
                       magnitudeRange = Range(6.0, 9.0),
                       depthRange = Range(0.0, 400.0),
                       distanceRange = None)

    print("\n")
    print(" ** Test ArcLink Request Builder **")
    print(" ** Station Based **")
    print("\n")
    arb = ArcLinkFDSNBuilder("http://www.moho.iag.usp.br", "seisrequest.iag.usp.br:18001:m.bianchi@iag.usp.br")
    re = arb.stationBased(t0 = UTCDateTime("2014-01-01"),
                          t1 = UTCDateTime("2015-01-01"),
                          targetSamplingRate = 20.0,
                          allowedLocGainList = ["H", "L"],
                          dataWindowRange = Range(-600.0, 600.0),
                          phasesOrPhaseGroupList = [ "pgroup" ],
                  
                          networkStationList = [ "BL.AQDB", "BL.PCMB" ],
                          stationRestrictionArea = AreaRange(-150.0, 60.0, 20.0, 55.0),
                  
                          eventRestrictionArea = None,
                          magnitudeRange = Range(7.0, 9.0),
                          depthRange = Range(0.0, 400.0),
                          distanceRange = Range(0, 90.0))

    print("\n")
    print(" ** Test Request Builder **")
    print(" ** Show Request **")
    print("\n")
    ArcLinkFDSNBuilder.show_request(rs, False)

    print("\n")
    print(" ** Test Request Builder **")
    print(" ** Show Request Compact **")
    print("\n")
    ArcLinkFDSNBuilder.show_request(rs)

    print("\n")
    print(" ** Test Fetching **")
    print(" ** ArcLink fetcher **")
    print("\n")
    fetcher = Sc3ArclinkFetcher("seisrequest.iag.usp.br:18001/m.bianchi@iag.usp.br",
                        allinone  = False,
                        merge     = True)

    dl = Downloader(basedir       = "./Test_Save/",
                      replacetree = True,
                      show_resume = False,
                      fetcher     = fetcher,
                      saverlist   = None)
    dl.work(rs)

    print("\n")
    print(" ** Test Saving (Sac + Qfile + MS modes 1,2,3) **")
    print(" ** ArcLink fetcher **")
    print("\n")

    assac = SacSaver()
    asq   = QSaver()

    asm1   = MSSaver(1) # Mode 1
    asm2   = MSSaver(2) # Mode 2
    asm3   = MSSaver(3) # Mode 3

    dl = Downloader(basedir       = "./Test_SaveE/",
                      replacetree = True,
                      show_resume = True,
                      fetcher     = fetcher,
                      saverlist   = [ assac, asq, asm1, asm2, asm3 ])
    dl.enableSaveRaw()
    dl.work(re)

    dl = Downloader(basedir = "./Test_SaveS/",
                      replacetree=True,
                      show_resume = True,
                      fetcher = fetcher,
                      saverlist = [ assac, asq, asm1, asm2, asm3 ])
    dl.work(rs)
