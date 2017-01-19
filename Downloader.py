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

import shutil, StringIO, os, sys
from seiscomp.arclink.manager import ArclinkManager, ArclinkError, arclink_status_string
from Savers import Saver

class BaseFetcher(object):
    pass

class FDSNFetcher(BaseFetcher):
    def __init__(self, hostorclient = "IRIS", allinone = False, merge = False):
        if isinstance(hostorclient, fdsn.Client):
            self._fc = hostorclient
            self._host = hostorclient.base_url
        elif isinstance(hostorclient, str):
            self._host = hostorclient
            self._fc = fdsn.Client(self._host)
        else:
            raise Exception("Invalid parameter hostorclient, type missmatch.")

        if allinone:
            print("This is only a compatibility parameter, will be ignored", file=sys.stderr)

        self._all_in_one = allinone
        self._merge = merge

    def work(self, key, items):
        stream = Stream()
        for (t0, t1, net, sta, channels, _, _, _) in items:
            for (loc, cha, _, _) in channels:
                try:
                    traces = self._fc.get_waveforms(net,sta,loc,cha,t0,t1)
                    if self._merge:
                        traces.merge(method=1, fill_value="interpolate")
                    stream += traces
                except fdsn.header.FDSNException,e:
                    print("%s.%s.%s.%s %s" % (net,sta,loc,cha,str(e)), file=sys.stderr)
                    continue
                except AttributeError,e:
                    # This is a fix for the client
                    print("%s.%s.%s.%s %s" % (net,sta,loc,cha,str(e)), file=sys.stderr)
                    continue
        return stream

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

        except ArclinkError,e:
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

class Downloader(object):
    ''' Downloader: 
        
        Should receive a folder (basedir), a flag to allow replace a tree (replacetree) a flag 
        to indicate that it should show a resume of operations (show_resume) and finally, should
        be given a fetcher (fetcher) of type (BaseFetcher) and a list of extracters (saverlist) of
        type (Savers).'''
    def __init__(self, basedir,
                 replacetree = False,
                 show_resume = False,
                 fetcher = None,
                 saverlist = None):
        self._basedir = None
        self._show_resume = show_resume
        self._force = replacetree
        self.__saveraw= False

        if fetcher is not None and not isinstance(fetcher, BaseFetcher):
            raise Exception("Supplied fetcher %s is not a Fetcher." % str(fetcher))
        self._fetcher = fetcher

        if not isinstance(saverlist, list):
            self._savers = [ saverlist ]
        else:
            self._savers = saverlist

        for saver in self._savers:
            if saver is not None and not isinstance(saver, Saver):
                raise Exception("Extracter %s is not of type Saver Class" % str(saver))

        if not os.path.isdir(basedir):
            try:
                os.mkdir(basedir)
            except OSError,e:
                print("Failed to initialize the base folder (%s)." % (str(e)), file=sys.stderr)
        self._basedir = basedir

    def isgood(self):
        if self._basedir == None: return False
        return True

    def _resume(self, key, items):
        i=0
        for item in items:
            print("  %03d) T1-T2: %s-%s Total: %d s Station: %s.%s - %d Channels" % (i,
                                                                                       item[0],
                                                                                       item[1],
                                                                                       item[1]-item[0],
                                                                                       item[2],
                                                                                       item[3],
                                                                                       len(item[4])), file=sys.stderr)
            i+=1

    def _buildfolder(self, key):
        return os.path.join(self._basedir, key)

    def _makefolders(self, requests):
        created = 0
        removed = 0

        # Loop on request inside requests
        for key in requests.keys():
            if key == "STATUS": continue

            # Find the folder
            _cfolder = self._buildfolder(key)

            # Check if exists
            if os.path.isdir(_cfolder):
                if not self._force:
                    raise Exception("Folder %s already exists" % _cfolder)
                shutil.rmtree(_cfolder)
                removed += 1

            # Create it
            os.mkdir(_cfolder)
            created += 1

        return (removed, created)

    def enableSaveRaw(self):
        folder = os.path.join(self._basedir, "RAW")
        if not os.path.isdir(folder): os.mkdir(folder)
        self.__saveraw = True

    def work(self, requests):

        print("\n\nWorking on %d requests" % (len(requests) - 1), file=sys.stderr)

        # Ensure folders exists
        (removed, created) = self._makefolders(requests)
        print(" Removed %d folders and created %d folder in %s" % (removed, created, self._basedir), file=sys.stderr)

        # Work on request base
        for key in requests.keys():
            if key == "STATUS": continue

            request = requests[key]
            print("\n %s has %d events selected" % (key,len(request)), file=sys.stderr)

            # Check resume
            if self._show_resume:
                self._resume(key, request)

            # Fetch
            if self._fetcher:
                data = self._fetcher.work(key, request)
                if data == None or len(data) == 0:
                    print("  No Data for Node %s" % key, file=sys.stderr)

                if data and self.__saveraw:
                    data.write(os.path.join(self._basedir, "RAW", "%s.mseed" % key), "MSEED")

                if data and self._savers:
                    for extracter in self._savers:
                        if extracter is None: continue
                        result = extracter.work(self._buildfolder(key), key, request, data)
                        print("  Wrote (In:%d Assoc:%d nWin:%d rms:%d spike:%d 3c:%d) -- %d files" % result, file=sys.stderr)

if __name__ == "__main__":
    from Builders import BaseBuilder
    rq = BaseBuilder.load_request("xc-phase2")
    A = Sc3ArclinkFetcher("seisrequest.iag.usp.br:18001:m.bianchi@iag.usp.br", allinone = True, merge = False)
    
    for key in rq.keys():
        if key == "STATUS": continue
        request = rq[key]
        print (request[0])
        A.work(key, [ request[0] ])
