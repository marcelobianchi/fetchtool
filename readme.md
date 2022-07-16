# FetchTool

FetchTool is a multi modular mass downloader tool for seismological data. It is aware of FDSN and ArcLink services for metadata and data download and saves data as SAC, Qfiles or MSeed today.

## Installing

To install it we recomend you to use pip along with a checkout of this repository since this is still in an alpha stage. It is being constant updated, to easy the constant upgrade you can install a link to your repo using 'pip install . -e', something like:

```
$ git clone https://github.com/marcelobianchi/fetchtool.git
$ cd fetchtool/
$ pip install . -e --user
```

> Don't forget to change pip to pip3 if needed. This code is not expected to work with python 2 anymore.

Also, if you want to fetch an update, just move to were you cloned it the first time and run git pull again. Something like:

```
$ cd fetchetool/
$ git pull
```

# Overview

This is a modular package. You have 3 different types of modules to build your program. The first module is the request builder, that can build request lines combining event and station in two different ways: by station and by events. The second module is the Downloader that reads the request data from the request builder and download waveform data from server. Finally, you have the savers that are given to the downloaders. Each saver is responsible to save the files in a different format and place. You have the SAC saver, a Qfile saver and in the future, a MS saver.

Also, it should be easy in the future to complement already mass downloaded data with new acquired data.

## Builders

There are two different builders, the FDSNBuilder that is a complete FDSN builder and an ArcLinkFDSNBuilder, that reads events from FDSN and metadata from ArcLink server.

## Downloader

The Downloader main module just download the data, but it needs the fetchers. There are available two different fetchers, a FDSN fetcher and an SCc3ArclinkFetcher that get data from an ArcLink server.

The Downloader module can be instructed to save the RAW data obtained from the server. Once data is downloaded from server it is checked and passed to each of the supplied saver module.

## Savers

Saver modules are modules that enforces constrains in the data and save it to certain pre-estabilished directory structure. Main savers implemented are the SacSaver, QSaver and a MSSaver. When saving the data they first fill in the headers with the maximum amount of available data. Also QSaver, prepares a SeismiHandler station file to easy the import of data into the program. MSSaver don't save any station or event information associated with the data.

# About the Author

The whole package was written by Marcelo Bianchi while working at the University of São Paulo. It uses the ObsPy package and ArcLink clients supplied by the SeisComP3 system. It is tested with ObsPy 1.0.

For comments and question please address directly Marcelo Bianchi (<m.bianchi@iag.usp.br>).

This is free software, under the GNU-GPL 3 license.

_Last Updated 24/March/2016_

