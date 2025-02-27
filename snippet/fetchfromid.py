from fetchtool.BaseBuilder import Range, USP_DATA, USP_EVENTS, USP_DATA_INTERNAL
from fetchtool.Builders import FDSNBuilder

from fetchtool.Savers import SacSaver
from fetchtool.Downloader import Downloader, FDSNFetcher

import sys

if __name__ == "__main__":
	EID = sys.argv[1]

	rb = FDSNBuilder(USP_EVENTS)
	rq = rb.eventidBased(EID, Range(-60, 360), stationDistanceRange = Range(0.0, 10.0))

	sv = SacSaver()
	ft = FDSNFetcher(USP_DATA_INTERNAL, allinone = True, merge = False)
	dl = Downloader("./", True, False, ft, [sv])

	done = []
	with open("snuffler.stations", "w") as fio:
		for l in rq:
			if l == 'STATUS' : continue
			lines = rq[l]
			for t0,t1,N,S,CLS,Sa,Ea,Pa in lines:
				if Sa.stationId in done:
					continue
				print(Sa.stationId + ".", Sa.latitude, Sa.longitude, Sa.elevation, 0.0, "-n/a-", file = fio)
				done.append(Sa.stationId)

	done = []
	with open("snuffler.events","w") as fio:
		for l in rq:
			if l == 'STATUS' : continue
			lines = rq[l]
			for t0,t1,N,S,CLS,Sa,Ea,Pa in lines:
				if Ea.eventId in done:
					continue
				print(f'name = {Ea.eventId}', file = fio)
				print(f'time = {Ea.time.strftime("%Y-%m-%d %H:%M:%S.%f")}', file = fio)
				print(f'latitude = {Ea.latitude}', file = fio)
				print(f'longitude = {Ea.longitude}', file = fio)
				print(f'magnitude = {Ea.magnitude}', file = fio)
				print(f'catalog = LOCAL', file = fio)
				done.append(Ea.eventId)

	dl.work(rq)
