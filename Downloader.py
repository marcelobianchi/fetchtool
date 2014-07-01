from obspy import fdsn
import shutil, StringIO, os, sys
from seiscomp.arclink.manager import ArclinkManager, ArclinkError
from Savers import Saver
from obspy.core import read as oREAD, Stream

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
            print >>sys.stderr,"This is only a compatibility parameter, will be ignored"

        self._all_in_one = allinone
        self._merge = merge

    def work(self, key, items):
        stream = Stream()
        for (t0, t1, net, sta, channels, sa, ev, evp) in items:
            for (loc, cha, a, d) in channels:
                try:
                    traces = self._fc.get_waveforms(net,sta,loc,cha,t0,t1)
                    if self._merge:
                        traces.merge(method=1, fill_value="interpolate")
                    stream += traces
                except fdsn.header.FDSNException,e:
                    print "%s.%s.%s.%s %s" % (net,sta,loc,cha,str(e))
                    continue
                except AttributeError,e:
                    # This is a fix for the client
                    print "%s.%s.%s.%s %s" % (net,sta,loc,cha,str(e))
                    continue
        return stream

class Sc3ArclinkFetcher(BaseFetcher):
    def __init__(self, hostorclient = "seisrequest.iag.usp.br:18001/anonymous@unreal.plc", allinone = False, merge = False):
        if isinstance(hostorclient, str):
            (self._host, self._user) = hostorclient.split("/")
            self._am = ArclinkManager(self._host, self._user)
        elif isinstance(hostorclient, ArclinkManager):
            self._am = hostorclient
            self._host = self._am._ArclinkManager__myaddr
            self._user = self._am._ArclinkManager__default_user
        else:
            raise Exception("Invalid object hostorclient -- type missmatch")

        self._all_in_one = allinone
        self._merge = merge

    def _logquery(self, status):
        message = "  Request issues:"

        for vol in status.volume:
            for line in vol.line:
                if line.status != 3:
                    if message:
                        print >>sys.stderr,message
                        message = None
                    print >>sys.stderr,"   erro on line %s reason %s" % (line.content,arclink_status_string(line.status))

    def _download(self, request, stream):
        datastream = StringIO.StringIO()

        # Submit the request
        request.submit(self._host, self._user)

        # Download the request
        try:
            request.wait()
            self._logquery(request.status())
            request.download_data(datastream, block = True, purge = True)

            if datastream.pos == 0:
                print >>sys.stderr,"  No Data Returned for Requests:"
                for line in request.content:
                    print >>sys.stderr,"  %s.%s.%s.%s %s %s" % (line.net, line.sta, line.loc, line.cha, line.start_time, line.end_time)
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
            print "  Arclink Message: %s" % str(e)

        return self._am.new_request("WAVEFORM")

    def work(self, key, items):
        stream = Stream()

        # Prepare the initial request
        rq = self._am.new_request("WAVEFORM")

        if self._merge: self._all_in_one = False

        # Go on
        for (t0, t1, net, sta, channels, sa, ev, evp) in items:

            for (location, channel, azimuth, dip) in channels:
                rq.add(t0.datetime, t1.datetime, net, sta, channel, location)

            if not self._all_in_one:
                rq = self._download(rq, stream)

        if self._all_in_one:
            if len(rq.content) > 400:
                raise Exception("Attention, you are trying to ask for more than 400 lines of request !!!!")

            self._download(rq, stream)

        return stream 

class Downloader(object):
    def __init__(self, basedir,
                 replacetree = False,
                 resume = False,
                 fetcher = None,
                 extracter = None):
        self._basedir = None
        self._show_resume = resume
        self._force = replacetree

        if fetcher is not None and not isinstance(fetcher, BaseFetcher):
            raise Exception("Supplied fetcher %s is not a Fetcher." % str(fetcher))
        self._fetcher = fetcher

        if not isinstance(extracter, list):
            self._extracter = [ extracter ]
        else:
            self._extracter = extracter

        for saver in self._extracter:
            if saver is not None and not isinstance(saver, Saver):
                raise Exception("Extracter %s is not of type Saver Class" % str(saver))

        if not os.path.isdir(basedir):
            try:
                os.mkdir(basedir)
            except OSError,e:
                print >>sys.stderr,"Failed to initialize the base folder (%s)." % (str(e))
        self._basedir = basedir

    def isgood(self):
        if self._basedir == None: return False
        return True

    def _resume(self, key, items):
        i=0
        for item in items:
            print >>sys.stderr,"  %03d) T1-T2: %s-%s Total: %d s Station: %s.%s - %d Channels" % (i,
                                                                                       item[0],
                                                                                       item[1],
                                                                                       item[1]-item[0],
                                                                                       item[2],
                                                                                       item[3],
                                                                                       len(item[4]))
            i+=1

    def _buildfolder(self, key):
        return os.path.join(self._basedir, key)

    def _makefolders(self, requests):
        created = 0
        removed = 0

        # Loop on request inside requests
        for key in requests.keys():

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

    def work(self, requests, extracterparams):

        print >>sys.stderr,"\n\nWorking on %d requests" % len(requests)

        # Ensure folders exists
        (removed, created) = self._makefolders(requests)
        print >>sys.stderr," Removed %d folders and created %d folder in %s" % (removed, created, self._basedir)

        # Work on request base
        for key in requests.keys():
            request = requests[key]
            print >>sys.stderr,"\n %s has %d events selected" % (key,len(request))

            # Check resume
            if self._show_resume:
                self._resume(key, request)

            # Fetch
            if self._fetcher:
                data = self._fetcher.work(key, request)
                if data == None or len(data) == 0:
                    print >>sys.stderr,"  No Data for Node %s" % key

                if data and self._extracter:
                    for extracter in self._extracter:
                        if extracter is None: continue
                        result = extracter.work(self._buildfolder(key), key, request, data, extracterparams)
                        print >>sys.stderr,"  Wrote (In:%d Assoc:%d nWin:%d rms:%d 3c:%d) -- %d files" % result 
