import re,os,sys

class ShStation(object):
    ''' Class to handle with STATINF.DAT files'''
    def __init__(self, filename, append = True):
        self._file = filename
        self._rawdata = []
        self._keys = {}
        self._modified = None

        if os.path.isfile(filename):
            if append:
                self._import()

    def _parseline(self,line,linec):
        items = {}
        stcode = None

        lloop = line
        ma = re.match("(?P<stcode>\w+)\s+((?P<name>\w+)[:]\s*(?P<val>[0-9.,+-]+))", lloop)

        if ma:
            stcode = ma.group("stcode")

            while ma != None:
                if ma.group("name") not in [ "lat", "lon", "elevation", "xrel", "yrel", "name", "array"]:
                    raise Exception("Invalid key name %s" % ma.group("name"))
                items[ma.group("name")] = ma.group("val") 
                lloop = lloop[ma.end():]
                ma = re.match("(\s+(?P<name>\w+)[:]\s*(?P<val>[a-zA-Z0-9.,+-_]+))", lloop)

            # Save the file linenumber
            items["linec"] = linec
            items["stcode"] = stcode

        return (stcode, items)

    def _import(self):
        stfile = open(self._file, "r")

        linec = 0
        for line in stfile:
            sline = line.strip()
            if len(sline) == 0 or sline[0] == '!':
                self._rawdata.append(line)
            else:
                (key, items) = self._parseline(line, linec)
                if key is None:
                    raise Exception("Invalid line %d: '%s'" % (linec,line))
                if key in self._keys:
                    raise Exception("Station is dupplicated !")
                self._keys[key.upper()] = linec
                self._rawdata.append(items)
            linec += 1
        stfile.close()
        self._modified = False

    def save(self, filename = None):
        if not filename and not self._modified: return
        print >>sys.stderr,"Saving file ... "
        if filename:
            iofile = open(filename,"w")
        else:
            if os.path.isfile(self._file):
                os.rename(self._file, self._file + "~")
            iofile = open(self._file, "w")

        for line in self._rawdata:
            if isinstance(line, str):
                iofile.write(line)
            elif isinstance(line,dict):
                iofile.write("%-5s " % line["stcode"])
                for k in [ "lat", "lon", "elevation", "array", "xrel", "yrel", "name"]:
                    if k in line:
                        iofile.write("%s:%s " % (k,line[k]))
                iofile.write("\n")

        if not filename:
            self._modified = False

        iofile.close()

    def add(self, stcode, lat, lon, elevation, array = None, xrel = None, yrel = None, name = None):
        if stcode in self._keys:
            raise Exception("Station is already defined")
        items = {}

        items["lat"] = lat
        items["lon"] = lon
        items["elevation"] = elevation

        # Optional items
        if array:
            items["array"] = array
        if xrel:
            items["xrel"] = xrel
        if yrel:
            items["yrel"] = yrel
        if name:
            items["name"] = name

        items["stcode"] = stcode
        items["linec"] = len(self._rawdata)
        self._keys[stcode] = len(self._rawdata)
        self._rawdata.append(items)
        self._modified = True

    def has(self, stcode):
        return stcode in self._keys

    def get(self,stcode):
        if stcode.upper() in self._keys:
            return self._rawdata[self._keys[stcode.upper()]]
        return None

