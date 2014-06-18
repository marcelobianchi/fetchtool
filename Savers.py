from SeismicHandlerStation import ShStation
import os, sys
from obspy.core import AttribDict

# Defaults Time Constants
SECOND = 1
MINUTE = 60*SECOND
HOUR = 60*MINUTE
DAY = 24*HOUR
WEEK = 7 * DAY
MONTH = 30 * DAY
HEAR = 365 * DAY

class Saver(object):
    def __init__(self, debug = False):
        self._debug = debug

    def _cleanids(self, reason, stream, ids):
        '''Remove from the stream the traces that have no evid tag
           or that have its _f_evid attribute listed in the ids
           variable'''

        trs = []

        if self._debug: print >>sys.stderr,"CleanIds called from %s" % reason
        while len(stream):
            trace = stream.pop()
            # Checks
            if not hasattr(trace, '_f_evid'):
                if self._debug: print >>sys.stderr," cleaning trace (%s/%s)" % (trace.id,trace.stats.starttime)
                continue
            if trace._f_evid in ids:
                if self._debug: print >>sys.stderr," cleaning trace (%s/%s) EvId: %s" % (trace.id,trace.stats.starttime, trace._f_evid)
                continue
            trs.append(trace)
        if self._debug: print >>sys.stderr,"Clean Done."

        stream.extend(trs)

    def _associate(self, stream, request):
        ''' Associate each trace to a single event in the list.'''

        extra = []

        # Tag traces that correspond to an phase arrival (associate event - trace)
        linecount = 0
        for (t0, t1, net, sta, channels, sa, ev, evp) in request:
            for trace in stream:
                ts = trace.stats.starttime
                te = trace.stats.endtime
                stid = "%s.%s" % (trace.stats.network,trace.stats.station)
                if evp['time'] > ts and evp['time'] <= te and stid == sa['stationId']:
                    if hasattr(trace, '_f_evid'):
                        if self._debug: print >>sys.stderr,"Duplicating trace %s evid=%s in favor of evid=%s" % (trace.id, trace._f_evid, ev['eventId'])
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

        # Collect number of channels per ids
        i = 0
        for trace in stream:
            try:
                id = ids[trace._f_evid]
            except KeyError:
                id= []
                ids[trace._f_evid] = id
            id.append((trace.stats.channel, i, trace.stats.npts))
            i += 1

        # Collect ids to remove
        ids_to_remove = [ ]
        for id in ids:
            zs = self._collect(ids[id], ["Z"])
            ns = self._collect(ids[id], ["N", "1"])
            es = self._collect(ids[id], ["E", "2"])

            zs = self._makechoice(stream, zs)
            ns = self._makechoice(stream, ns)
            es = self._makechoice(stream, es)

            if len(zs) != 1 or len(ns) != 1 or len(es) != 1:
                ids_to_remove.append(id)

        # Clean-up
        self._cleanids("3C-", stream, ids_to_remove)

    def _check_minimun_window_size(self, stream, request, parameters):
        if not parameters.tw.prephasevalue and not parameters.tw.postphasevalue: return

        # Search for bad IDS
        for trace in stream:
            ts = trace.stats.starttime
            te = trace.stats.endtime
            pt = request[trace._f_linecount][7]['time']
            if (pt - ts) < parameters.tw.prephasevalue:
                del trace._f_evid
                continue
            if (te - pt) < parameters.tw.postphasevalue:
                del trace._f_evid
                continue

        # Clean up bad IDS
        self._cleanids("Minimun Window", stream, [])

    def _ensure_minimun_rms(self, stream, request, parameters):
        pass

    def _fix_event_headers(self, stream, request):
        raise Exception("Base Class Saver -- Not implemented")

    def _fix_station_headers(self, stream, request):
        raise Exception("Base Class Saver -- Not implemented")

    def _extract(self, folder, key, request, stream):
        raise Exception("Base Class Saver -- Not implemented")

    def getDefaultParameters(self):
         # TimeWindow
        parameters = AttribDict({})
        parameters.tw = AttribDict({})
        parameters.tw.prephasevalue = 1 * MINUTE
        parameters.tw.postphasevalue = 10 * MINUTE

        # RMS
        parameters.rms = AttribDict({})
        parameters.rms.rmsratio = 2.0
        parameters.rms.prephasevalue = 5
        parameters.rms.postphasevalue = 25
        parameters.rms.freqmin = None
        parameters.rms.freqmax = None

        return parameters

    def work(self, folder, key, request, stream, parameters):
        if not os.path.isdir(folder): raise Exception("Invalid folder")

        n_initial = len(stream)

        ## Process Stream Traces
        self._associate(stream, request)
        n_associate = len(stream)

        # Garantee that all traces has at least ...
        self._check_minimun_window_size(stream, request, parameters)
        n_window = len(stream)

        # Garantee that the phase has an RMS of at least ...
        self._ensure_minimun_rms(stream, request, parameters)
        n_rms = len(stream)

        ## THIS IS THE LAST
        # Ensure that every event has 3 components only !
        self._ensure_tree_components(stream, request)
        n_tree = len(stream)

        # Fill the headers available for SAC file format
        self._fix_station_headers(stream, request)
        self._fix_event_headers(stream, request)

        ## Extract
        written = self._extract(folder, key, request, stream)

        return (n_initial, n_associate, n_window, n_rms, n_tree, written)

class QSaver(Saver):
    def __init__(self, debug = False):
        Saver.__init__(self, debug)

    def _getfilename(self, key):
        return "%s" % (key.upper())

    def _extract(self, folder, key, request, stream):
        filename = self._getfilename(key)
        if self._debug: print >>sys.stderr,"  Wrote %s" % os.path.join(folder, filename)
#        stream.write("debug", format="MSEED")
        stream.sort(['EVIDMINE','channel'])
        stream.write(folder, format="Q")

        # Handle the STATINF.DAT SH file
        try:
            statinf = ShStation(os.path.join(folder,"..","STATINF.DAT"))
            for item in request:
                stp = item[5]
                stcode = stp.stationId.split(".")[1]
                if not statinf.has(stcode):
                    statinf.add(stcode, stp.latitude, stp.longitude, stp.elevation)
            statinf.save()
            del statinf
        except Exception,e:
            print str(e)

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

            trace.stats.sh['P-ONSET'] = phase.time

    def _fix_station_headers(self, stream, request):
        ''' Seismic Handler does not store info 
            about stations in the file ! '''
        pass

class SacSaver(Saver):
    def __init__(self, debug = False):
        Saver.__init__(self, debug)

    def _getfilename(self,trace):
        return "%s_%s.sac" % (trace.id, trace._f_evid)

    def _extract(self, folder, key, request, stream):
        nw = 0
        for trace in stream:
            filename = self._getfilename(trace)
            if self._debug: print >>sys.stderr,"  Wrote %s" % os.path.join(folder, filename)
            trace.write(os.path.join(folder, filename), format="SAC")
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
#             trace.stats.sac['mag'] = None

            # Selection phase
            trace.stats.sac['a'] = phase.time - trace.stats.starttime
            trace.stats.sac['ka'] = phase.phase

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
                    trace.stats.sac['cmpaz'] = cha[2]
                    trace.stats.sac['cmpinc'] = cha[3] + 90.0
                    break
            else:
                print >>sys.stderr,"Cannot decide on channel %s orientation" % trace.id
                print >>sys.stderr,"File %s will not have orientation set" % (self._getfilename(trace))
