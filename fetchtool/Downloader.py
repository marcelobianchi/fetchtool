'''Downloader orchestration Module and Fetchers for the fetchtool package.

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
from fetchtool.Savers import Saver
from socket import timeout

class BaseFetcher(object):
    '''A Fetcher super class
    '''
    def work(self, key, items):
        '''The main fetcher worker
        
        When implementing a new fetcher you are responsible to implement this.
        
        Parameters
        ----------
        key : str
            The key of items
        items: list
            A lista of request lines corresponding to the request key received
        
        Return
        ------
        obspy.Stream
            An obspy stream with the fetched data
        
        '''
        raise Exception("Not implemented")

class FDSNFetcher(BaseFetcher):
    '''A fetcher based on FDSN server
    
    Parameters
    ----------
    hostorclient : str, default "IRIS"
        The URL of the host to connect to. See the warning above if you need to supply a user and password.
    allinone : bool, default False
        Indicate that the whole request should be request as one request.
    merge, bool, default False
        Indicate that before returning the data to the Downloader the fetcher should merge the received traces.
    
    Warning
    -------
        If you want to obtain restricted data you will need to indicate a user and a password on the hostorclient.
        Another option is that you instantiate yourself an obspyclient and supply it to the FDSNFetcher. The current
        synthax to indicate the **username** and **password** is using
        
        hostorclient = "http://YOUR_FDSN_HOST;YOUR_USER;YOUR_PASSWORD"
        
    '''

    def __init__(self, hostorclient = "IRIS", allinone = False, merge = False):
        if isinstance(hostorclient, fdsn.Client):
            print("Using supplied client at: {}".format(hostorclient.base_url))
            self._fc = hostorclient
        elif isinstance(hostorclient, str):
            host = hostorclient
            user = None
            password = None
            
            if ";" in hostorclient:
                host,user,password = hostorclient.split(";", 2)
            
            self._fc = fdsn.Client(host, user = user, password = password)
        else:
            raise Exception("Invalid parameter hostorclient, type missmatch.")

        if self._fc.user is None:
            print("Connected to {} without authentication.".format(self._fc.base_url))
        else:
            print("Connected to {} authenticated as user={}".format(self._fc.base_url, self._fc.user))
        
        
        if allinone:
            print("Since all in one is given will BULK all requests into one", file=sys.stderr)

        self._all_in_one = allinone
        self._merge = merge

    def work(self, key, items):
        '''This is the main method of the FDSNFetcher.
        
        This is used by the Downloadeder. You should not call it yourself.
        
        '''
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
                print("%s.%s %s" % (bulk[0][0],bulk[0][1],str(e)), file=sys.stderr)
            except AttributeError as e:
                # This is a fix for the client
                print("%s.%s %s" % (bulk[0][0],bulk[0][1],str(e)), file=sys.stderr)
            except timeout as e:
                print("%s.%s %s" % (bulk[0][0],bulk[0][1],str(e)), file=sys.stderr)
            except:
                print("%s.%s Unknown error" % (bulk[0][0],bulk[0][1]), file=sys.stderr)
        
        return stream

class RAWFetcher(BaseFetcher):
    '''A fetcher that reads files from the RAW folder.
    
    This fetcher can be used to (re)load the data that was once saved in a prior
    run that you need to re-process without downloading it again. It is also very
    usefull for debuging.
    
    Parameters
    ----------
    basedir : str
        The base folder where the files are, or, that has a RAW folder in it with files.

    '''
    
    def __init__(self, basedir):
        self._basedir = basedir

    def work(self, key, items):
        filename = os.path.join(self._basedir, "RAW", "%s.mseed" % key)
        altfilename = os.path.join(self._basedir, "%s.mseed" % key)
        
        if os.path.isfile(filename):
            return oREAD(filename)
        elif os.path.isfile(altfilename):
            return oREAD(altfilename)
        
        return Stream()

class Downloader(object):
    '''The Downloader
        
        This is the main downloader module. It will save all the requested data into the basedir folder,
        if requested it can overwrite already downloaded data (replacetree). One downloader work with one
        fetcher and can use many different Savers.
        
        Parameters
        ----------
        basedir : str
            A base folder path to save files.
        replacetree : bool, default False
            A flag to indicate that the folder should be replaced (overwrite files)
        show_resume: bool, default False
            A flag to indicate that a print-out of the request should be shown
        fetcher : BaseFetcher instance
            A fetcher instance.
        saverlist: Saver list
            A list of Saver type objects.
        '''
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

    ''' Public Methods '''
    
    def isgood(self):
        '''A check method for the Downloader'''
        
        if self._basedir == None:
            return False
        
        if os.path.isdir(self._basedir) and self._force is False:
            return False
        
        return True

    def enableSaveRaw(self):
        '''Request to save RAW data before QC
        '''
        folder = os.path.join(self._basedir, "RAW")
        if not os.path.isdir(folder): os.mkdir(folder)
        self.__saveraw = True

    def work(self, request):
        '''Main download method.
        
        This method will receive a request object, download the data and save it to the base folder
        
        Parameters
        ----------
        request : request
            A request to work on
        
        Return
        ------
        None
        '''
        
        print("\n\nWorking on %d request entries" % (len(request) - 1), file=sys.stderr)

        # Ensure folders exists
        (removed, created) = self._makefolders(request)
        print(" Removed %d folders and created %d folder in %s" % (removed, created, self._basedir), file=sys.stderr)

        # Work on request base
        for key in request.keys():
            if key == "STATUS": continue

            reqitem = request[key]
            print("\n %s has %d lines selected" % (key, len(reqitem)), file=sys.stderr)

            # Check resume
            if self._show_resume:
                self._resume(key, reqitem)

            # Fetch
            if self._fetcher:
                data = self._fetcher.work(key, reqitem)
                if data == None or len(data) == 0:
                    print("  No Data for Node %s" % key, file=sys.stderr)

                if data and self.__saveraw:
                    data.write(os.path.join(self._basedir, "RAW", "%s.mseed" % key), "MSEED")

                if data and self._savers:
                    for extracter in self._savers:
                        if extracter is None: continue
                        result = extracter.work(self._buildfolder(key), key, reqitem, data)
                        print("  Wrote (In:%d Assoc:%d nWin:%d rms:%d spike:%d 3c:%d nHead:%d) -- %d files" % result, file=sys.stderr)

if __name__ == "__main__":
    pass
