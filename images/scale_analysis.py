#!/usr/local/bin/python

# This script performs the scale analysis between IMERG/TMPA and MRMS, 
# generating the files needed by the plotting scripts.


import numpy as np
import gzip
import h5py
import os
import pickle
from glob import glob
from calendar import monthrange
from datetime import datetime, timedelta
from pyhdf.SD import SD, SDC


datapath = '/home/jackson/Work/Data/Scale Analysis/'
mrmspath = '/media/jackson/Vault/MRMS/Level3/'
imergpath = '/media/jackson/Vault/IMERG-F/'
tmpapath = '/media/jackson/Vault/TRMM/3B42/'

year0, year1 = 2014, 2015
month0, month1 = 6, 12
nday = (datetime(year1, month1, monthrange(year1, month1)[1]) - 
        datetime(year0, month0, 1)).days + 1
daterange = '%4d%02d-%4d%02d' % (year0, month0, year1, month1)

times = (1, 3, 6, 12, 24)
sizes1 = list(range(1, 26))
sizes2 = list(range(1, 11))
ensemble = 100

# first, define some convenience function

def getMonths(year, year0, month0, year1, month1):
    if   year0 == year1: return range(month0, month1 + 1)
    elif year  == year0: return range(month0, 13)
    elif year  == year1: return range(1, month1 + 1)
    else               : return range(1, 13)


# Part 1: We extract MRMS rainfall within 30°N−41.5°N, 93.5°W−83.5°W and write
# them to file. Note that these files retain the MRMS dimensions of (hour, lat, 
# lon), with latitude running from north to south.


lat0 = int((54.995 - 41.495) / 0.01)
lat1 = int((54.995 - 30.005) / 0.01 + 1)
lon0 = int((-93.495 + 129.995) / 0.01)
lon1 = int((-83.505 + 129.995) / 0.01 + 1)
nlat, nlon = lat1 - lat0, lon1 - lon0


for year in range(year0, year1 + 1):
    for month in getMonths(year, year0, month0, year1, month1):
        for day in range(1, monthrange(year, month)[1] + 1):

            outfile1 = '%sUSEast/p_%4d%02d%02d.npy.gz' % (mrmspath, year, month, day)
            outfile2 = '%sUSEast/r_%4d%02d%02d.npy.gz' % (mrmspath, year, month, day)
            outfile3 = '%sUSEast/s_%4d%02d%02d.npy.gz' % (mrmspath, year, month, day)

            if (not os.path.exists(outfile1) or 
                not os.path.exists(outfile2) or
                not os.path.exists(outfile3)):

                P = np.zeros([48, nlat, nlon], dtype = np.float32)
                R = np.zeros([48, nlat, nlon], dtype = np.float32)
                S1 = np.zeros([48, nlat, nlon], dtype = np.float32)
                S2 = np.zeros([48, nlat, nlon], dtype = np.float32)

                for tt in range(48):

                    T = datetime(year, month, day, tt // 2, (tt % 2) * 30) + \
                            timedelta(minutes = 30)
                    yy, mm, dd, HH, MM = T.year, T.month, T.day, T.hour, T.minute

                    inpath = '%s%4d/%02d/' % (mrmspath, yy, mm)
                    ts = '%4d%02d%02d.%02d%02d00' % (yy, mm, dd, HH, MM)
                    infile1 = '30MGCP.%s.asc.gz' % (ts,)
                    infile2 = '30MRQI.%s.asc.gz' % (ts,)
                    infile3 = '30MTYPE03.%s.asc.gz' % (ts,)
                    infile4 = '30MTYPE04.%s.asc.gz' % (ts,)

                    try:
                        with gzip.open('%s%s' % (inpath, infile1), 'rb') as f:
                            data = f.readlines()[6:]
                        P[tt] = np.array([[float(ii) for ii in jj.split()[lon0 : lon1]] 
                                           for jj in data[lat0 : lat1]])
                        with gzip.open('%s%s' % (inpath, infile2), 'rb') as f:
                            data = f.readlines()[6:]
                        R[tt] = np.array([[float(ii) for ii in jj.split()[lon0 : lon1]] 
                                           for jj in data[lat0 : lat1]])
                        with gzip.open('%s%s' % (inpath, infile3), 'rb') as f:
                            data = f.readlines()[6:]
                        S1[tt] = np.array([[float(ii) for ii in jj.split()[lon0 : lon1]] 
                                            for jj in data[lat0 : lat1]])
                        with gzip.open('%s%s' % (inpath, infile4), 'rb') as f:
                            data = f.readlines()[6:]
                        S2[tt] = np.array([[float(ii) for ii in jj.split()[lon0 : lon1]] 
                                            for jj in data[lat0 : lat1]])
                    except FileNotFoundError:
                        pass

                # mask values using NaN
                P[P < 0] = np.nan
                R[R < 0] = np.nan
                S = S1 + S2
                S[S < 0] = np.nan

                with gzip.open(outfile1, 'wb') as f:
                    np.save(f, P)
                with gzip.open(outfile2, 'wb') as f:
                    np.save(f, R)
                with gzip.open(outfile3, 'wb') as f:
                    np.save(f, S)


# Part 2: To get the mean rain rates for IMERG, TMPA and MRMS, we begin by 
# selecting the location of the boxes, with one box for each size and ensemble. 
# We record the top-left corner of the boxes, making sure that the entire box
# remains within our region. This is done separately for IMERG and TMPA. We 
# then calculate the mean IMERG and MRMS rain rates as well as the mean TMPA 
# and MRMS rain rates.


# generate the coordinates for the IMERG boxes

inpath = '%s%4d/%02d/%02d/' % (imergpath, 2014, 4, 1)
files = sorted(glob('%s3B-HHR*' % inpath))
nt = len(files)
with h5py.File(files[0]) as f:
    lats = f['Grid/lat'][:]
    lons = f['Grid/lon'][:]
    fillvalue = f['Grid/precipitationCal'].attrs['_FillValue']
lat0 = np.where(lats ==  30.05)[0][0]
lat1 = np.where(lats ==  41.45)[0][0] + 1
lon0 = np.where(lons == -93.45)[0][0]
lon1 = np.where(lons == -83.55)[0][0] + 1

file_lat = '%slat_IMERG_%s.npy.gz' % (datapath, daterange)
file_lon = '%slon_IMERG_%s.npy.gz' % (datapath, daterange)

if os.path.exists(file_lat) and os.path.exists(file_lon):

    with gzip.open(file_lat, 'rb') as f:
        imerg_lat = np.load(f)
    with gzip.open(file_lon, 'rb') as f:
        imerg_lon = np.load(f)

else:
    
    imerg_lat = np.zeros([len(sizes1), ensemble], dtype = np.int32)
    imerg_lon = np.zeros([len(sizes1), ensemble], dtype = np.int32)

    for ss, size in enumerate(sizes1):
        lat = np.random.randint(lat0, lat1 - size, ensemble)
        lon = np.random.randint(lon0, lon1 - size, ensemble)
        imerg_lat[ss] = np.array(lat)
        imerg_lon[ss] = np.array(lon)
        
    with gzip.open(file_lat, 'wb') as f:
        np.save(f, imerg_lat)
    with gzip.open(file_lon, 'wb') as f:
        np.save(f, imerg_lon)

# calculate IMERG for the boxes

Pimerg = {}
for size in sizes1:
    for en in range(ensemble):
        Pimerg[0.5, size, en] = np.full([nday, nt], np.nan, dtype = np.float32)

dor = 0   # day of record
for year in range(year0, year1 + 1):
    for month in getMonths(year, year0, month0, year1, month1):
        for day in range(1, monthrange(year, month)[1] + 1):

            inpath = '%s%4d/%02d/%02d/' % (imergpath, year, month, day)
            files = sorted(glob('%s3B-HHR*' % inpath))

            for tt in range(nt):
                with h5py.File(files[tt], 'r') as f:
                    for ss, size in enumerate(sizes1):
                        for en in range(ensemble):
                            lat = imerg_lat[ss, en]
                            lon = imerg_lon[ss, en]
                            P = f['Grid/precipitationCal'][lon : lon + size, 
                                                           lat : lat + size]
                            mask = P == fillvalue
                            if not np.all(mask):
                                Pimerg[0.5, size, en][dor, tt] = np.mean(P[~mask])

            dor += 1

file_imerg = '%sIMERG_%s.gz' % (datapath, daterange)
with gzip.open(file_imerg, 'wb') as f:
    pickle.dump(Pimerg, f)

# calculate MRMS for the IMERG boxes

Pmrms1 = {}
for size in sizes1:
    for en in range(ensemble):
        Pmrms1[0.5, size, en] = np.full([nday, 48], np.nan, dtype = np.float32)

dor = 0
for year in range(year0, year1 + 1):
    for month in getMonths(year, year0, month0, year1, month1):
        for day in range(1, monthrange(year, month)[1] + 1):

            infile1 = '%sUSEast/p_%4d%02d%02d.npy.gz' % (mrmspath, year, month, day)
            infile2 = '%sUSEast/r_%4d%02d%02d.npy.gz' % (mrmspath, year, month, day)
            infile3 = '%sUSEast/s_%4d%02d%02d.npy.gz' % (mrmspath, year, month, day)

            with gzip.open(infile1) as f:
                P = np.load(f)[:, ::-1]
            with gzip.open(infile2) as f:
                R = np.load(f)[:, ::-1]
            with gzip.open(infile3) as f:
                S = np.load(f)[:, ::-1]

            for tt in range(48):
                for ss, size in enumerate(sizes1):
                    for en in range(ensemble):
                        lat = (imerg_lat[ss, en] - lat0) * 10
                        lon = (imerg_lon[ss, en] - lon0) * 10

                        p = P[tt, lat : lat + size * 10, lon : lon + size * 10]
                        r = R[tt, lat : lat + size * 10, lon : lon + size * 10]
                        s = S[tt, lat : lat + size * 10, lon : lon + size * 10]
                        mask = np.isnan(p) + np.isnan(r) + np.isnan(s)
                        if not np.all(mask):
                            Pmrms1[0.5, size, en][dor, tt] = np.mean(p[~mask])

            dor += 1

file_mrms1 = '%sMRMS1_%s.gz' % (datapath, daterange)
with gzip.open(file_mrms1, 'wb') as f:
    pickle.dump(Pmrms1, f)

# generate the coordinates for the TMPA boxes

lons = np.arange(-179.875, 180, 0.25)
lats = np.arange(-49.875, 50, 0.25)

lat0 = np.where(lats ==  30.125)[0][0]
lat1 = np.where(lats ==  41.375)[0][0] + 1
lon0 = np.where(lons == -93.375)[0][0]
lon1 = np.where(lons == -83.625)[0][0] + 1

file_lat = '%slat_TMPA_%s.npy.gz' % (datapath, daterange)
file_lon = '%slon_TMPA_%s.npy.gz' % (datapath, daterange)

if os.path.exists(file_lat) and os.path.exists(file_lon):

    tmpa_lat = np.load(file_lat)
    tmpa_lon = np.load(file_lon)

else:
    
    tmpa_lat = np.zeros([len(sizes2), ensemble], dtype = np.int32)
    tmpa_lon = np.zeros([len(sizes2), ensemble], dtype = np.int32)

    for ss, size in enumerate(sizes2):
        lat = np.random.randint(lat0, lat1 - size, ensemble)
        lon = np.random.randint(lon0, lon1 - size, ensemble)

        tmpa_lat[ss] = np.array(lat)
        tmpa_lon[ss] = np.array(lon)
        
    with open(file_lat, 'wb') as f:
        np.save(f, tmpa_lat)
    with open(file_lon, 'wb') as f:
        np.save(f, tmpa_lon)

# calculate TMPA for the boxes

Ptmpa = {}
for size in sizes2:
    for en in range(ensemble):
        Ptmpa[3, size, en] = np.full([nday, 8], np.nan, dtype = np.float32)

dor = 0   # day of record
for yy1 in range(year0, year1 + 1):
    for mm1 in getMonths(yy1, year0, month0, year1, month1):
        for dd1 in range(1, monthrange(yy1, mm1)[1] + 1):
            for hh, hh1 in enumerate(range(0, 24, 3)):

                # next time step
                ts = datetime(yy1, mm1, dd1, hh1) + timedelta(hours = 3)
                yy2, mm2, dd2, hh2 = ts.year, ts.month, ts.day, ts.hour

                # unzip the files if they are not available
                inpath1 = '%s%4d/%02d/%02d/' % (tmpapath, yy1, mm1, dd1)
                infile1 = '%s3B42.%4d%02d%02d.%02d.7.HDF' % (inpath1, yy1, mm1, 
                                                             dd1, hh1)
                if not os.path.exists(infile1):
                    with open(infile1, 'wb') as f:
                        with gzip.open('%s.gz' % infile1) as g:
                            f.write(g.read())

                inpath2 = '%s%4d/%02d/%02d/' % (tmpapath, yy2, mm2, dd2)
                infile2 = '%s3B42.%4d%02d%02d.%02d.7.HDF' % (inpath2, yy2, mm2, 
                                                             dd2, hh2)
                if not os.path.exists(infile2):
                    with open(infile2, 'wb') as f:
                        with gzip.open('%s.gz' % infile2) as g:
                            f.write(g.read())

                f1 = SD(infile1, SDC.READ)
                P1 = f1.select('precipitation').get()
                f2 = SD(infile2, SDC.READ)
                P2 = f2.select('precipitation').get()

                for ss, size in enumerate(sizes2):
                    for en in range(ensemble):

                        lon = tmpa_lon[ss, en]
                        lat = tmpa_lat[ss, en]

                        P = (P1[lon : lon + size, lat : lat + size] + 
                             P2[lon : lon + size, lat : lat + size]) / 2
                        mask = P < -1
                        if not np.all(mask):
                            Ptmpa[3, size, en][dor, hh] = np.mean(P[~mask])

                f1.end()
                f2.end()

                os.remove(infile1)

            dor += 1

os.remove(infile2)

file_tmpa = '%sTMPA_%s.gz' % (datapath, daterange)
with gzip.open(file_tmpa, 'wb') as f:
    pickle.dump(Ptmpa, f)

# calculate MRMS for the TMPA boxes

Pmrms2 = {}
for size in sizes2:
    for en in range(ensemble):
        Pmrms2[0.5, size, en] = np.full([nday, 48], np.nan, dtype = np.float32)

dor = 0
for year in range(year0, year1 + 1):
    for month in getMonths(year, year0, month0, year1, month1):
        for day in range(1, monthrange(year, month)[1] + 1):

            infile1 = '%sUSEast/p_%4d%02d%02d.npy.gz' % (mrmspath, year, month, day)
            infile2 = '%sUSEast/r_%4d%02d%02d.npy.gz' % (mrmspath, year, month, day)
            infile3 = '%sUSEast/s_%4d%02d%02d.npy.gz' % (mrmspath, year, month, day)

            with gzip.open(infile1) as f:
                P = np.load(f)[:, ::-1]
            with gzip.open(infile2) as f:
                R = np.load(f)[:, ::-1]
            with gzip.open(infile3) as f:
                S = np.load(f)[:, ::-1]

            for tt in range(48):
                for ss, size in enumerate(sizes2):
                    for en in range(ensemble):
                        lat = (tmpa_lat[ss, en] - lat0) * 25
                        lon = (tmpa_lon[ss, en] - lon0) * 25

                        p = P[tt, lat : lat + size * 25, lon : lon + size * 25]
                        r = R[tt, lat : lat + size * 25, lon : lon + size * 25]
                        s = S[tt, lat : lat + size * 25, lon : lon + size * 25]
                        mask = np.isnan(p) + np.isnan(r) + np.isnan(s)
                        if not np.all(mask):
                            Pmrms2[0.5, size, en][dor, tt] = np.mean(p[~mask])

            dor += 1

file_mrms2 = '%sMRMS2_%s.gz' % (datapath, daterange)
with gzip.open(file_mrms2, 'wb') as f:
    pickle.dump(Pmrms2, f)
