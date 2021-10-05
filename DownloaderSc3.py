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

from __future__ import division, print_function

from obspy.clients import fdsn
from obspy import read as oREAD, Stream

try:
    import StringIO
except:
    from io import StringIO

import shutil, os, sys
from seiscomp.arclink.manager import ArclinkManager, ArclinkError, arclink_status_string
from Savers import Saver

class Sc3ArclinkFetcher(BaseFetcher):
    def __init__(self, hostorclient = "seisrequest.iag.usp.br:18001:anonymous@unreal.plc", allinone = False, merge = False):
        if isinstance(hostorclient, str):
            (host, port, self._user) = hostorclient.replace("/",":").split(":")
            self._host = "%s:%s" % (host,port)
            self._am = ArclinkManager(self._host, self._user)
        elif isinstance(hostorclient, ArclinkManager):
            self._am = hostorclient
            self._host = self._am._ArclinkManager__myaddr
            self._user = self._am._ArclinkManager__default_user
        else:
            raise Exception("Invalid object hostorclient -- type missmatch")

        self._all_in_one = allinone
        self._merge = merge
        self.__SSLpasswordDict = None

    def __load_passwords(self, extrafile  = None):
        if self.__SSLpasswordDict is not None:
            return
        
        self.__SSLpasswordDict = {}
        
        for SSLpasswordFile in [ os.path.join(os.environ['HOME'],".dcidpasswords.txt"), "dcidpasswords.txt", extrafile]: 
            if SSLpasswordFile and os.path.exists(SSLpasswordFile):
                fd = open(SSLpasswordFile)
                line = fd.readline()
                while line:
                    line = line.strip()
                    if line and line[0] != "#":
                        try:
                            (dcid, password) = line.split()
                            if dcid in self.__SSLpasswordDict:
                                print("Overwritting %s datacenter default password with password from file %s" % (dcid, SSLpasswordFile))
                            self.__SSLpasswordDict[dcid] = password
                        except ValueError:
                            print(SSLpasswordFile + " invalid line: " + line, file=sys.stderr)
                    line = fd.readline()
                fd.close()

    def _logquery(self, status):
        message = "  Request issues:"

        for vol in status.volume:
            for line in vol.line:
                if line.status != 3:
                    if message:
                        print(message,file=sys.stderr)
                        message = None
                    print("   erro on line %s reason %s" % (line.content,arclink_status_string(line.status)), file=sys.stderr)

    def _download(self, request, stream):
        datastream = StringIO.StringIO()

        # Submit the request
        request.submit(self._host, self._user)
        if request.error:
            print("  Error while submiting Arclink Request: %s" % request.error, file=sys.stderr)
            return self._am.new_request("WAVEFORM")
        
        # Download the request
        try:
            request.wait()
            self._logquery(request.status())
            
            dcid = None
            for v in request.status().volume:
                if v.encrypted:
                    self.__load_passwords()
                    if dcid and dcid != v.dcid:
                        print("Oops, volumes from different dcids ... %s - not supported" % v.dcid)
                    dcid = v.dcid
            
            if self.__SSLpasswordDict is not None and dcid and dcid in self.__SSLpasswordDict:
                request.download_data(datastream, block = True, purge = True, password = self.__SSLpasswordDict[dcid])
            else:
                request.download_data(datastream, block = True, purge = True)

            if datastream.pos == 0:
                print("  No Data Returned for Requests:", file=sys.stderr)
                for line in request.content:
                    print("  %s.%s.%s.%s %s %s" % (line.net, line.sta, line.loc, line.cha, line.start_time, line.end_time), file=sys.stderr)
                return self._am.new_request("WAVEFORM")

            datastream.seek(0)
            traces = oREAD(datastream)
            del datastream
    
            if self._merge:
                traces.merge(method=1, fill_value="interpolate")
    
            # Append Streams
            stream += traces

        except ArclinkError as e:
            self._logquery(request.status())
            request.purge()
            print("  Arclink Message: %s" % str(e))

        return self._am.new_request("WAVEFORM")

    def work(self, key, items):
        stream = Stream()

        # Prepare the initial request
        rq = self._am.new_request("WAVEFORM")

        if self._merge: self._all_in_one = False

        # Go on
        for (t0, t1, net, sta, channels, _, _, _) in items:

            for (location, channel, _, _) in channels:
                rq.add(t0.datetime, t1.datetime, net, sta, channel, location)

            if not self._all_in_one or len(rq.content) > 400:
                if len(rq.content) > 400:
                    print("Downloading first %d lines of request" % len(rq.content), file=sys.stderr)
                rq = self._download(rq, stream)

        if self._all_in_one:
            if len(rq.content) > 400:
                print("Attention, you are trying to ask for more than 400 lines of request !!!!", file=sys.stderr)

            self._download(rq, stream)

        return stream

