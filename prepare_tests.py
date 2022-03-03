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

from Builders import FDSNBuilder
from BaseBuilder import Range, AreaRange

if __name__ == "__main__":
    
    Bb = FDSNBuilder("IRIS", "USP")

    req = Bb.eventBased(t0 = UTCDateTime("2017-01-01"), t1 = UTCDateTime("2018-01-01"), 
                  targetSamplingRate = 100, allowedLocGainList = [ "H" ],
                  dataWindowRange = Range(5,20), phasesOrPhaseGroupList =  "ttp", eventRestrictionArea = AreaRange.WORLD(), magnitudeRange = Range(8,10), 
                  depthRange = Range.ALLDEPTHS(), networkStationList = [ "BL.AQDB", "BL.PCMB" ])

    Bb.save_request("test_request.rq", req, True)
