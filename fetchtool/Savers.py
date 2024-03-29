'''A set of Savers to be used with fetchtool

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

from __future__ import print_function, division

from fetchtool.SeismicHandlerStation import ShStation
from obspy.core import AttribDict, Stream
from obspy.io.sac import SACTrace
from obspy.geodetics import gps2dist_azimuth, kilometers2degrees
import os, sys
import numpy as np 
import warnings

# Defaults Time Constants
SECOND = 1
''' A second'''
MINUTE = 60  * SECOND
''' Minute constant in seconds'''
HOUR   = 60  * MINUTE
''' Hour constant in seconds'''
DAY    = 24  * HOUR
''' Day constant in seconds'''
WEEK   = 7   * DAY
''' Week constant in seconds'''
MONTH  = 30  * DAY
''' Month constant in seconds'''
YEAR   = 365 * DAY
''' Year constant in seconds'''

'''
Base Saver
'''
class Saver(object):
    '''The base superclass for Savers
    '''
    def __init__(self, debug = False):
        self._debug = debug
        self.parameters = self.__initParameters()
        '''Set of current SAVER parameter'''

    def _cleanids(self, reason, stream, ids):
        '''Remove from the stream the traces that have no evid tag
           or that have its _f_evid attribute listed in the ids
           variable'''

        trs = []

        if self._debug: print("CleanIds called from %s" % reason, file=sys.stderr)
        while len(stream):
            trace = stream.pop()
            # Checks
            if not hasattr(trace, '_f_evid'):
                if self._debug: print(" cleaning trace (%s/%s)" % (trace.id,trace.stats.starttime), file=sys.stderr)
                continue
            if trace._f_evid in ids:
                if self._debug: print(" cleaning trace (%s/%s) EvId: %s" % (trace.id,trace.stats.starttime, trace._f_evid), file=sys.stderr)
                continue
            trs.append(trace)
        if self._debug: print("Clean Done.", file=sys.stderr)

        stream.extend(trs)

    def _associate(self, stream, request):
        ''' Associate each trace to a single event in the list.'''

        extra = []

        # Tag traces that correspond to an phase arrival (associate event - trace)
        linecount = 0
        for (_, _, _, _, _, sa, ev, evp) in request:
            for trace in stream:
                ts = trace.stats.starttime
                te = trace.stats.endtime
                stid = "%s.%s" % (trace.stats.network,trace.stats.station)
                if evp['time'] > ts and evp['time'] <= te and stid == sa['stationId']:
                    if hasattr(trace, '_f_evid'):
                        if self._debug: print("Duplicating trace %s evid=%s in favor of evid=%s" % (trace.id, trace._f_evid, ev['eventId']), file=sys.stderr)
                        trace = trace.copy()
                        extra.append(trace)
                    trace._f_evid = "%s_%s" % (ev['eventId'], sa['stationId'])
                    trace.stats.EVIDMINE = ev['eventId']
                    trace._f_linecount = linecount
            linecount += 1

        stream.extend(extra)

        # Clean-up
        self._cleanids("Associate", stream, [])

    def _collect(self, items, pattern):
        gooditems = []
        for (cha, n, npts) in items:
            if cha[2] in pattern:
                gooditems.append((cha, n, npts))
        return gooditems

    def _makechoice(self, stream, items):
        if len(items) == 1: return items

        nptsmax = None
        gooditem = (None, None, None)
        for item in items:
            if nptsmax is None or nptsmax <= item[2]:
                nptsmax = item[2]
                gooditem = item

        for item in items:
            if item == gooditem: continue
            del stream[item[1]]._f_evid

        if gooditem == (None,None,None):
            return []
        else:
            return [ gooditem ]

    def _reqmode(self, request, Net, Sta):
        '''List the channels orientation that are listed in request for network and station.
        
        Parameters
        ----------
        request : tuple
            A request tuple
        Net : str
            The network code string
        Sta : str
            The station code string
        
        Return
        ------
        str
            A string with the orientations listed in request for net and sta.
        '''
        for (_, _, N, S, Cs, _, _, _) in request:
            if N == Net and S == Sta:
                return "".join([ cha[-1] for (loc, cha, _, _) in Cs ])
        return ""

    def _ensure_all_components(self, stream, request):
        ids = {}
        req_ids = {}
        if self.parameters.nocompc: return

        # Collect number of channels per ids
        i = 0
        for trace in stream:
            try:
                evid = ids[trace._f_evid]
            except KeyError:
                evid= []
                ids[trace._f_evid] = evid
                req_ids[trace._f_evid] = self._reqmode(request, trace.stats.network, trace.stats.station)
            # BUG DEVERIA INCLUIR A LOCATION
            evid.append((trace.stats.channel, i, trace.stats.npts))
            i += 1

        # Collect ids to remove
        ids_to_remove = [ ]
        for evid in ids:
            zs = self._collect(ids[evid], ["Z"])
            ns = self._collect(ids[evid], ["N", "1"])
            es = self._collect(ids[evid], ["E", "2"])

            zs = self._makechoice(stream, zs)
            ns = self._makechoice(stream, ns)
            es = self._makechoice(stream, es)

            if "Z" in req_ids[evid] and len(zs) != 1:
                ids_to_remove.append(evid)
                continue

            if ("N" in req_ids[evid] or "1" in req_ids[evid]) and len(ns) !=1:
                ids_to_remove.append(evid)
                continue

            if ("E" in req_ids[evid] or "2" in req_ids[evid]) and len(es) != 1:
                ids_to_remove.append(evid)
                continue

        ids_to_remove = list(set(ids_to_remove))

        # Clean-up
        self._cleanids("3C-", stream, ids_to_remove)

    def _check_minimun_window_size(self, stream, request):
        if not self.parameters.tw.prephasevalue and not self.parameters.tw.postphasevalue: return

        # Search for bad IDS
        for trace in stream:
            ts = trace.stats.starttime
            te = trace.stats.endtime
            pt = request[trace._f_linecount][7]['time']
            if pt - ts < self.parameters.tw.prephasevalue:
                del trace._f_evid
                continue
            if te - pt < self.parameters.tw.postphasevalue:
                del trace._f_evid
                continue

        # Clean up bad IDS
        self._cleanids("Minimun Window", stream, [])

    def _ensure_no_spike(self, stream, request):
        if self.parameters.spike.enabled == False: return

        for trace in stream:
            pt = request[trace._f_linecount][7]['time']

            ta = trace.slice(pt - self.parameters.spike.window, pt)
            ta.detrend('linear')

            tb = trace.slice(pt, pt + self.parameters.spike.window)
            tb.detrend('linear')

            if abs((ta.max() - ta.min())/2.) * self.parameters.spike.ratio > abs((tb.max() - tb.min())/2.):
                del trace._f_evid

        # Clean up bad IDS
        self._cleanids("Spike Clean", stream, [])

    def _ensure_minimun_rms(self, stream, request):
        if self.parameters.rms.enabled == False: return
        
        for trace in stream:
            pt = request[trace._f_linecount][7]['time']

            ta = trace.slice(pt-self.parameters.rms.w, pt)
            ta.detrend('linear')

            tb = trace.slice(pt, pt+self.parameters.rms.w)
            tb.detrend('linear')

            rmsa = np.mean(ta.data * ta.data)
            rmsb = np.mean(tb.data * tb.data)

            if rmsa * self.parameters.rms.ratio >= rmsb:
                del trace._f_evid

        # Clean up bad IDS
        self._cleanids("RMS Check", stream, [])

    @staticmethod
    def _decimate_resolver(from_sps, to_sps):
        table = {}
        done = [ 0 ]
        work = [ to_sps ]

        # Build a reduction table
        ##
        while max(done) <= from_sps:
            base = work.pop(0)
            if base in done: continue
            table[base] = dict([ (i*base,[i]) for i in [2,3,5,7]])
            work.extend(list(table[base].keys()))
            work.sort()
            done.append(base)

        coefs = []
        ctable = table[to_sps]
        del table[to_sps]

        # Reduce
        ##
        while table:
            if from_sps in ctable:
                return ctable[from_sps] # Found it
            
            for k,kl in list(ctable.items()):
                if k not in table: continue
                for i,il in list(table[k].items()):
                    if i in ctable: continue
                    ctable[i] = il + kl
                del table[k]

        return None

    def _resample_data(self, stream, request):
        if not self.parameters.resample.enabled: return
        
        for trace in stream:
            if trace.stats.sampling_rate == self.parameters.resample.sps: continue
            chain = self._decimate_resolver(trace.stats.sampling_rate, self.parameters.resample.sps)
            
            if self._debug:
                print('Decimation from {} to {} using {}'.format(trace.stats.sampling_rate, self.parameters.resample.sps, (chain if chain is not None else "INTERPOLATION")), file = sys.stderr)

            if chain is None:
                if trace.stats.sampling_rate >= self.parameters.resample.sps:
                    maxfreq =  (self.parameters.resample.sps / 2.) * 0.95
                    print('Filtering prior to interpolation f < {} to match sps = {}'.format(maxfreq, self.parameters.resample.sps), file = sys.stderr)
                    trace.filter('lowpass_cheby_2', freq = maxfreq)
                
                trace.interpolate(sampling_rate = self.parameters.resample.sps)
                continue
            
            for link in chain:
                trace.decimate(link)
            
        return

    def _fix_event_headers(self, stream, request):
        raise Exception("Base Class Saver -- Not implemented")

    def _fix_station_headers(self, stream, request):
        raise Exception("Base Class Saver -- Not implemented")

    def _extract(self, folder, key, request, stream):
        '''This is the main class of a new Saver.
        
        When writting a new saver you are responsible to implement this
        method. This method will receive the destination folder, the key
        of the request corresponding to the Stream, the request element
        and the stream with data.
        
        You should write it to the folder.
        
        Parameters:
        
        folder: str
            Where you should write the data to
        
        key : str
            The key of the request being saved
        
        request : request
            The whole request, request[key] is the one being saved
        
        stream : obspy.Stream
            The data to be saved.
        
        Return
        ------
        int
            The number of files written
        '''
        raise Exception("Base Class Saver -- Not implemented")

    def enableTimeWindowCheck(self, prephase, postphase):
        '''Enable a time window check for the downloaded data.
        
        When this option is enabled the saver will not save traces that 
        does not have at leasr prephase and postphase seconds before and
        after the reference phase.
        
        Parameters
        ----------
        prephase: float
            Amount of seconds before the reference phase (value is considered to be absolute)
        postphase: float
            Amount of seconds after the reference phase (value is considered to be absolute)
        '''
        
        self.parameters.tw.prephasevalue = abs(prephase)
        self.parameters.tw.postphasevalue = abs(postphase)

    def enablermscheck(self, window, ratio):
        '''When enabled will check SNR (by the RMS) of data before save
        
        Any traces with SNR ratio < than ratio will not be saved.
        
        Parameters
        ----------
        window : float
            The window size in seconds before and after the reference phase to check
        ratio:float
            The minimum ratio to save data windows
        '''
        self.parameters.rms.enabled = True
        self.parameters.rms.w = window
        self.parameters.rms.ratio = ratio

    def disablecompcheck(self):
        '''Disabling the check of components availability.
        
        Normaly requests that does not have all components are discarded.
        If disabled, data for stations missing requested components will be 
        accepted to be further processed.  
        '''
        self.parameters.nocompc = True

    def disable3ccheck(self):
        '''See the disablecompcheck method
        '''
        warnings.warn("This is depracated; use disablecompcheck", DeprecationWarning, stacklevel = 2)
        
        self.disablecompcheck()

    def enablespikecheck(self, window, ratio):
        '''Enable large amplitude check before the reference phase
        
        When enabled it will reject requests where any of traces has 
        a large amplitude before the reference fase. Diferently than
        the rms check this will use the (max-min)/2 ABS amplitude in the
        window before and after the reference phase to discard bad data.  
        
        Parameters
        ----------
        window : float
            A window in seconds to check for large amplitudes.
        
        ratio : float
            A minimum ratio to save the data.
        '''
        
        self.parameters.spike.enabled = True
        self.parameters.spike.window  = window
        self.parameters.spike.ratio   = ratio

    def disableresample(self):
        '''Disable resampling after it has been enabled with enableresample
        '''
        self.parameters.resample.enabled = False
        self.parameters.resample.sps = None

    def enableresample(self, sps):
        '''Enable data resample to sps before saving
        
        Parameters
        ----------
        sps : float
            A new sampling rate value to ressample data to.
        '''
        self.parameters.resample.enabled = True
        
        if sps < 1:
            print('SPS value = {} < 1, assuming dt was given!'.format(sps), file = sys.stderr)
            sps = int(1/sps)
        
        self.parameters.resample.sps = sps
        if self._debug: print('Enabling decimate/interpolate using an SPS = {}.'.format(sps), file=sys.stderr)
    
    def __initParameters(self):
        parameters = AttribDict({})

        # TimeWindow
        parameters.tw = AttribDict({})
        parameters.tw.prephasevalue = None
        parameters.tw.postphasevalue = None

        # RMS
        parameters.rms = AttribDict({})
        parameters.rms.enabled = False
        parameters.rms.ratio   = None
        parameters.rms.w       = None

        # 3c check
        parameters.nocompc = False

        # Spike
        parameters.spike = AttribDict({})
        parameters.spike.enabled = False
        parameters.spike.window  = None
        parameters.spike.ratio   = None

        # Ressampler
        parameters.resample = AttribDict({})
        parameters.resample.enabled = False
        parameters.resample.sps = None
        
        return parameters

    def work(self, folder, key, request, stream):
        '''The main working method.
        
        You should not call this directly. It is used by the downloader.
        '''
        if not os.path.isdir(folder): raise Exception("Invalid folder")

        stream = stream.copy()

        n_initial = len(stream)

        ## Process Stream Traces
        self._associate(stream, request)
        n_associate = len(stream)

        # Garantee that all traces has at least ...
        self._check_minimun_window_size(stream, request)
        n_window = len(stream)

        # Garantee that the phase has an RMS of at least ...
        self._ensure_minimun_rms(stream, request)
        n_rms = len(stream)

        self._ensure_no_spike(stream, request)
        n_spike = len(stream)
        
        ## THIS IS THE LAST
        # Ensure that every event has all the requested components !
        self._ensure_all_components(stream, request)
        n_tree = len(stream)

        # Fill the headers available for SAC file format
        self._fix_event_headers(stream, request)
        self._fix_station_headers(stream, request)
        n_head = len(stream)


        ## Resample & Extract
        self._resample_data(stream, request)
        written = self._extract(folder, key, request, stream)

        del stream

        return (n_initial, n_associate, n_window, n_rms, n_spike, n_tree, n_head, written)


'''
Savers
'''
class QSaver(Saver):
    '''The Q-File format Saver
    
    Parameters
    ----------
    debug : bool, default False
        Indicate that it should show debug messages while working.
    
    usenet_inname : bool, default False
        This indicates that the network code should be made part of the station code
    '''
    def __init__(self, debug = False, usenet_inname = False):
        Saver.__init__(self, debug)
        self.__usenet_inname = usenet_inname

    def _extract(self, folder, key, request, stream):
        if len(stream) == 0: return 0
        base, target = os.path.split(folder)
        if self._debug: print("  Wrote %s" % os.path.join(folder, "%s.Q??" % target), file=sys.stderr)

        stream.sort(['EVIDMINE','network', 'station', 'channel'])
        stream.write(os.path.join(folder, target), format="Q")

        ##
        # Handle the STATINF.DAT SH file
        ##
        try:
            statinf = ShStation(os.path.join(base,"STATINF.DAT"))
            for item in request:
                stp = item[5]
                if self.__usenet_inname:
                    stcode = "".join(stp.stationId.split("."))
                else:
                    stcode = stp.stationId.split(".")[1]
                if not statinf.has(stcode):
                    statinf.add(stcode, stp.latitude, stp.longitude, stp.elevation)
            statinf.save()
            del statinf
        except Exception as e:
            print(str(e))

        return 1

    def _fix_event_headers(self, stream, request):
        for trace in stream:
            evp = request[trace._f_linecount][6]
            phase = request[trace._f_linecount][7]
            if not hasattr(trace.stats,"sh"):
                trace.stats.sh = AttribDict({})

            trace.stats.sh['ORIGIN'] = evp.time
            trace.stats.sh['LAT'] = evp.latitude
            trace.stats.sh['LON'] = evp.longitude
            trace.stats.sh['DEPTH'] = evp.depth
            trace.stats.sh['COMMENT'] = evp.eventId

            try:
                if evp.magnitude is not None:
                        trace.stats.sh['MAGNITUDE'] = evp.magnitude
            except KeyError:
                print("No magnitude value (%s) set.", file=sys.stderr)

            trace.stats.sh['P-ONSET'] = phase.time

    def _fix_station_headers(self, stream, request):
        stream_head = Stream()
        stream_aux = Stream()
        aux = []

        for trace in stream:
            stp = request[trace._f_linecount][5]
            chalist = request[trace._f_linecount][4]
            if not hasattr(trace.stats,"sh"):
                trace.stats.sh = AttribDict({})

            # Station code
            if self.__usenet_inname:
                trace.stats.sh['STATION'] = "".join(stp.stationId.split("."))
            else:
                trace.stats.sh['STATION'] = stp.stationId.split(".")[1]

            # Orientation
            for cha in chalist:
                sid = "%s.%s.%s" % (stp.stationId, cha[0], cha[1])
                if sid == trace.id:
                    trace.stats.sh['DCVREG']  = cha[2]
                    trace.stats.sh['DCVINCI'] = cha[3]
                    stream_aux += trace
                    aux.append('OK')
                    break
            else:
                print("Cannot decide on channel %s orientation" % trace.id, file=sys.stderr)
                aux.append('NOT OK')
            
            if len(aux)==3:
                if 'NOT OK' not in aux:
                    stream_head += stream_aux
                    stream_aux = Stream()
                    aux = []
                else:
                    print("Will not use %s ... moving on\n" % (stp.stationId), file=sys.stderr)
                    stream_aux = Stream()
                    aux = []

        return stream_head


class SacSaver(Saver):
    '''The Sac-file format saver
    
    This is the SACS file format saver. It will load all information
    that it can inside the SAC header trace.
    
    It will also compute distance, gcarc, baz and az.
    
    Parameters
    ----------
    debug : bool, default False
        Debug while working to save data
    '''
    def __init__(self, debug = False):
        Saver.__init__(self, debug)

    def _getfilename(self, trace):
        try:
            filename = "%s.%s.%s.%s.sac" % (trace.stats.network, trace.stats.station, trace.stats.sac['kevnm'], trace.stats.channel)
        except:
            filename = "%s_%s.sac" % (trace.id, trace._f_evid)
        return filename

    def _extract(self, folder, key, request, stream):
        nw = 0
        for trace in stream:
            filename = self._getfilename(trace)
            if self._debug: print("  Wrote %s" % os.path.join(folder, filename), file=sys.stderr)
            # This is the way to go but does not work because of 
            # a bug in SAC module from ObsPy
            #trace.write(os.path.join(folder, filename), format="SAC")

            # HACK TO FIX OBSpy problems with SAC
            # http://lists.swapbytes.de/archives/obspy-users/2016-February/001991.html
            sac = SACTrace.from_obspy_trace(trace)
            for k in trace.stats.sac:
                setattr(sac, k, trace.stats.sac[k])
            sac.write(os.path.join(folder, filename))

            nw += 1
        return nw

    def _fix_event_headers(self, stream, request):
        for trace in stream:
            stp = request[trace._f_linecount][5]
            evp = request[trace._f_linecount][6]
            phase = request[trace._f_linecount][7]
            if not hasattr(trace.stats,"sac"):
                trace.stats.sac = AttribDict({})

            # Event parameters
            trace.stats.sac['o'] = evp.time - trace.stats.starttime
            trace.stats.sac['evla'] = evp.latitude
            trace.stats.sac['evlo'] = evp.longitude
            trace.stats.sac['evdp'] = evp.depth
            trace.stats.sac['kevnm'] = evp.eventId

            try:
                if evp.magnitude is not None:
                    trace.stats.sac['mag'] = evp.magnitude
            except KeyError as e:
                print("No magnitude value (%s)." % e, file=sys.stderr)

            # Distance
            distance, baz, az = gps2dist_azimuth(stp.latitude, stp.longitude, evp.latitude, evp.longitude)
            trace.stats.sac['dist'] = distance / 1000.0
            trace.stats.sac['gcarc'] = kilometers2degrees(distance / 1000.0)
            trace.stats.sac['baz'] = baz
            trace.stats.sac['az'] = az

            # Selection phase
            trace.stats.sac['a'] = phase.time - trace.stats.starttime
            trace.stats.sac['ka'] = phase.phase
            
            if phase.others is not None:
                for which,var in [ ('S','t0'), ('Pg','t1'),  ('Pn','t2'), ('Sg','t3'),  ('Sn','t4') ]:
                    for other in phase.others:
                        if other.phase != which: continue
                        trace.stats.sac[var] = other.time - trace.stats.starttime
                        trace.stats.sac['k' + var] = other.phase

    def _fix_station_headers(self, stream, request):
        for trace in stream:
            stp = request[trace._f_linecount][5]
            chalist = request[trace._f_linecount][4]

            if not hasattr(trace.stats,"sac"):
                trace.stats.sac = AttribDict({})

            # Station Coordinates
            trace.stats.sac['stla'] = stp.latitude
            trace.stats.sac['stlo'] = stp.longitude
            trace.stats.sac['stel'] = stp.elevation

            # Other settings
            trace.stats.sac['lcalda'] = True

            # Orientation
            for cha in chalist:
                sid = "%s.%s.%s" % (stp.stationId, cha[0], cha[1])
                if sid == trace.id:
                    trace.stats.sac['cmpaz']  = cha[2]
                    trace.stats.sac['cmpinc'] = cha[3] + 90.0
                    break
            else:
                print("Cannot decide on channel %s orientation" % trace.id, file=sys.stderr)
                print("File %s will not have orientation set" % (self._getfilename(trace)), file=sys.stderr)


class MSSaver(Saver):
    '''This class implements a miniSeed saver.

    Constructed filenames will be determined from the operaion mode.
    Possible values are:
    
     * 1 = File per key-channel (like sac) in key folder
     * 2 = File per key-station (like a sac with 3c) in key folder
     * 3 = File per key like Qfiles
    
    Parameters
    ----------
    mode : int
        The mode of operation
    debug : bool, default False
        Print debuging messages
    '''
    def __init__(self, mode, debug = False):
        Saver.__init__(self, debug)
        
        self._mode = mode
        
        if mode not in [1, 2, 3]:
            raise Exception("Invalid mode, choose between 1,2,3")
    
    def _fix_event_headers(self, stream, request):
        pass
    
    def _fix_station_headers(self, stream, request):
        pass
    
    def _getfilename(self, trace):
        if self._mode == 1:
            return "%s_%s.ms" % (trace.id, trace._f_evid)
        elif self._mode == 2:
            n, s, _, _ = trace.id.split(".")
            return "%s.%s_%s.ms" % (n, s, trace._f_evid)
        elif self._mode == 3:
            return None
    
    def _extract(self, folder, key, request, stream):
        nw = 0
        base, target = os.path.split(folder)
        if self._mode == 1:
            for trace in stream:
                filename = self._getfilename(trace)
                trace.write(os.path.join(folder, filename), format="MSEED")
                nw += 1

        elif self._mode == 2:
            streams = { }

            for trace in stream:
                filename = self._getfilename(trace)
                if filename not in streams:
                    streams[filename] = Stream()
                streams[filename] += trace

            for filename in streams:
                streams[filename].write(os.path.join(base, target, filename), format="MSEED")
                nw += 1

        elif self._mode == 3:
            stream.write(os.path.join(base, "%s.ms" % target), format="MSEED")
            nw = 1

        return nw
