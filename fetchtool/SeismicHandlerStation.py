'''A class to represent SeismicHandler station information

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
import re, os

class ShStation(object):
    '''Class to handle Seismic Handler station data inside a STATINF.DAT like file.
    
    Parameters
    ----------
    filename : str
        The filename of the file to handle.
    append : bool default True
        A flag to indicate if new stations should be appended to the file.
    '''
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
        '''Save the data to the original open file or to a new file
        
        Parameters
        ----------
        filename : str 
            The optional filename to save data to. If it is None, the
            filename indicated on class instance is used.
        '''
        if not filename and not self._modified: return
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
                iofile.write("%-7s " % line["stcode"])
                for k in [ "lat", "lon", "elevation", "array", "xrel", "yrel", "name"]:
                    if k in line:
                        iofile.write("%s:%s " % (k,line[k]))
                iofile.write("\n")

        if not filename:
            self._modified = False

        iofile.close()

    def add(self, stcode, lat, lon, elevation, array = None, xrel = None, yrel = None, name = None):
        '''Add a station record
        
        This method add a station record to the module. If it already in record an exception is raised.
        
        Parameters
        ----------
        stcode: str
            The station code
        lat : float
            The station latitude in degrees
        lon : float
            The station lonigtude in degrees
        elevation : float
            The station elevation in meters
        array : int, default None
            The number of the array
        xrel : float, default None
            The xoffset inside the array
        yrel : float, default None
            The yoffset inside the array
        name : str, default None
            The name/description of the station
        
        '''
        if stcode in self._keys:
            raise Exception("Station is already defined")
        items = {}

        items["lat"] = float(lat)
        items["lon"] = float(lon)
        items["elevation"] = float(elevation)

        # Optional items
        if array:
            items["array"] = array
        if xrel:
            items["xrel"] = float(xrel)
        if yrel:
            items["yrel"] = float(yrel)
        if name:
            items["name"] = str(name)

        items["stcode"] = stcode
        items["linec"] = len(self._rawdata)
        self._keys[stcode] = len(self._rawdata)
        self._rawdata.append(items)
        self._modified = True

    def has(self, stcode):
        '''Check that stcode is defined
        
        Parameters
        ----------
        stcode : str
            The station code to check
        
        Return
        ------
        bool
            True (exists) or False (don't exist)
        '''
        return stcode in self._keys

    def get(self,stcode):
        '''Returns a record from the class
        
        Parameters
        ----------
        stcode : str
            The code of the station to fetch
        
        Return
        ------
        dict
            A dictionary with station information
        '''
        if stcode.upper() in self._keys:
            return self._rawdata[self._keys[stcode.upper()]]
        
        return None

