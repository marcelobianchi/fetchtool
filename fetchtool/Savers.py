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

from __future__ import print_function, division

from fetchtool.SeismicHandlerStation import ShStation
from obspy.core import AttribDict, Stream
from obspy.io.sac import SACTrace
import os, sys
import numpy as np 

# Defaults Time Constants
SECOND = 1
MINUTE = 60  * SECOND
HOUR   = 60  * MINUTE
DAY    = 24  * HOUR
WEEK   = 7   * DAY
MONTH  = 30  * DAY
HEAR   = 365 * DAY

class Saver(object):
    def __init__(self, debug = False):
        self._debug = debug
        self.parameters = self.__initParameters()

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

    def _ensure_tree_components(self, stream, request):
        ids = {}
        if self.parameters.no3c: return

        # Collect number of channels per ids
        i = 0
        for trace in stream:
            try:
                evid = ids[trace._f_evid]
            except KeyError:
                evid= []
                ids[trace._f_evid] = evid
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

            if len(zs) != 1 or len(ns) != 1 or len(es) != 1:
                ids_to_remove.append(evid)

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

            if abs(ta.max()) * self.parameters.spike.ratio > abs(tb.max()):
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
                print(f'Decimation from {trace.stats.sampling_rate} to {self.parameters.resample.sps} using {chain if chain is not None else "INTERPOLATION"}', file = sys.stderr)
            
            if chain is None:
                if trace.stats.sampling_rate >= self.parameters.resample.sps:
                    maxfreq =  (self.parameters.resample.sps / 2.) * 0.95
                    print(f'Filtering prior to interpolation f < {maxfreq} to match sps = {self.parameters.resample.sps}', file = sys.stderr)
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
        raise Exception("Base Class Saver -- Not implemented")

    def enableTimeWindowCheck(self, prephase, postphase):
        self.parameters.tw.prephasevalue = abs(prephase)
        self.parameters.tw.postphasevalue = abs(postphase)

    def enablermscheck(self, window, ratio):
        self.parameters.rms.enabled = True
        self.parameters.rms.w = window
        self.parameters.rms.ratio = ratio

    def disable3ccheck(self):
        self.parameters.no3c = True

    def enablespikecheck(self, window, ratio):
        self.parameters.spike.enabled = True
        self.parameters.spike.window  = window
        self.parameters.spike.ratio   = ratio

    def disableresample(self):
        self.parameters.resample.enabled = False
        self.parameters.resample.sps = None

    def enableresample(self, sps):
        self.parameters.resample.enabled = True
        
        if sps < 1:
            print(f'SPS value = {sps} < 1, assuming dt was given!', file = sys.stderr)
            sps = int(1/sps)
        
        self.parameters.resample.sps = sps
        if self._debug: print(f'Enabling decimate/interpolate using an SPS = {sps}.', file=sys.stderr)
    
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
        parameters.no3c = False

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
        # Ensure that every event has 3 components only !
        self._ensure_tree_components(stream, request)
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

class QSaver(Saver):
    def __init__(self, debug = False, usenet_inname = False):
        Saver.__init__(self, debug)
        self.usenet_inname = usenet_inname

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
                if self.usenet_inname:
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
            if self.usenet_inname:
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
    def __init__(self, debug = False):
        Saver.__init__(self, debug)

    def _getfilename(self,trace):
        return "%s_%s.sac" % (trace.id, trace._f_evid)

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
    '''
    This class implements a miniSeed saver. It can
    save mseed in three different flavors indicated
    by the mode variable:
     * 1 : File per key-channel (like sac) in key folder
     * 2 : File per key-station (like a sac with 3c) in key folder
     * 3 : File per key like Qfiles
    '''
    def __init__(self, mode, debug = False):
        Saver.__init__(self, debug)
        self._mode = mode
        if mode not in [1, 2, 3]:
            raise Exception("Invalid mode, choose between 1,2,3 !")
    
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