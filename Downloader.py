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

import shutil, os, sys
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
            print("Since all in one is given will BULK all requests into one", file=sys.stderr)

        self._all_in_one = allinone
        self._merge = merge

    def work(self, key, items):
        stream = Stream()
        
        allbulk = []
        bulk = []
        
        for (t0, t1, net, sta, channels, _, _, _) in items:
            for (loc, cha, _, _) in channels:
                bulk.append((net,sta,loc,cha,t0,t1))
            
            if self._all_in_one:
                continue
            
            allbulk.append(bulk)
            bulk = []
        
        if len(bulk) > 0:
            allbulk.append(bulk)
        
        for bulk in allbulk:
            try:
                traces = self._fc.get_waveforms_bulk(bulk)
                if self._merge:
                    traces.merge(method=1, fill_value="interpolate")
                stream += traces
            except fdsn.header.FDSNException as e:
                print("%s.%s.%s.%s %s" % (net,sta,loc,cha,str(e)), file=sys.stderr)
            except AttributeError as e:
                # This is a fix for the client
                print("%s.%s.%s.%s %s" % (net,sta,loc,cha,str(e)), file=sys.stderr)
            except timeout as e:
                print("%s.%s.%s.%s %s" % (net,sta,loc,cha,str(e)), file=sys.stderr)
            except:
                print("%s.%s.%s.%s Unknown error" % (net,sta,loc,cha), file=sys.stderr)
        
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
            except OSError as e:
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
    pass
