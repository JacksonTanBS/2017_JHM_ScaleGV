#!/usr/local/bin/python

# This script provides the supporting functions to the plotting script.


import numpy as np
import matplotlib.pyplot as plt
import brewer2mpl
import os
import sys
import pickle
import gzip
from calendar import monthrange
from datetime import datetime

# setting up the parameters
year0, year1 = 2014, 2015
month0, month1 = 6, 12
nday = (datetime(year1, month1, monthrange(year1, month1)[1]) - 
        datetime(year0, month0, 1)).days + 1
times1 = (0.5, 1, 3, 6, 12, 24) # IMERG times (hours) to calculate
sizes1 = list(range(1, 26))     # IMERG sizes to calculate
times2 = (3, 6, 12, 24)         # TMPA times (hours) to calculate
sizes2 = list(range(1, 11))     # TMPA sizes to calculate
x1 = np.array(sizes1) * 0.1     # convert IMERG sizes to units of degree
x2 = np.array(sizes2) * 0.25    # convert TMPA sizes to units of degree
ensemble = 100                  # ensemble size

# plot configurations
cols1 = brewer2mpl.get_map('Set1', 'Qualitative', '7').hex_colors
del(cols1[5])
cols2 = cols1[2:]
scol = 3.503    # single column (89 mm)
dcol = 7.204    # double column (183 mm)
flpg = 9.724    # full page length (247 mm)
plt.rcParams['font.size'] = 9
plt.rcParams['legend.fontsize'] = 8
plt.rcParams['axes.titlesize'] = 'medium'
plt.rcParams['savefig.bbox'] = 'tight'
plt.rcParams['font.sans-serif'] = ['TeX Gyre Heros', 'Helvetica',
                                   'Bitstream Vera Sans']
plt.rcParams['pdf.fonttype'] = 42
dashes = (3, 3)    # dash spacing


def accum(data_base, res, base = 1):

    '''Function to calculate rate for longer periods from the base 
       rate (default: base = 1 h). Input units must be in mm / h and not 
       mm (i.e. not accumulation). Note: data_base must be in dimensions of 
       (day, base) and in units of mm / base.'''

    if (24 / res) % 1:

        print('Error: specified period not a factor of 24.')
        return False

    else:

        step = int(res / base)    # no. of steps to average over

        x = [np.nansum(data_base[:, ii * step : (ii + 1) * step], 1)
             for ii in range(24 // res)]

        return np.ma.array(x).T / step


def read_data(fixed_thresholds = False, season = 'all'):

    datapath = '/home/jackson/Work/Data/Scale Analysis/'
    daterange = '%4d%02d-%4d%02d' % (year0, month0, year1, month1)
    file_imerg = '%sIMERG_%s.gz' % (datapath, daterange)
    file_mrms1 = '%sMRMS1_%s.gz' % (datapath, daterange)
    file_tmpa  = '%sTMPA_%s.gz' % (datapath, daterange)
    file_mrms2 = '%sMRMS2_%s.gz' % (datapath, daterange)

    if fixed_thresholds:
        thresholds1, thresholds2 = {}, {}
        for time in times:
            for size in sizes1:
                thresholds1[time, size] = 0.1
            for size in sizes2:
                thresholds2[time, size] = 0.1
    else:
        try:
            thresholds1 = pickle.load(open('data_thresholds1.p', 'rb'))
            thresholds2 = pickle.load(open('data_thresholds2.p', 'rb'))
        except:
            print('Warning: thresholds files missing. Run threshold to generate.')
            thresholds1 = None
            thresholds2 = None

    Pimerg = np.load(gzip.open(file_imerg, 'rb'))
    Pmrms1 = np.load(gzip.open(file_mrms1, 'rb'))
    Ptmpa  = np.load(gzip.open(file_tmpa, 'rb'))
    Pmrms2 = np.load(gzip.open(file_mrms2, 'rb'))

    # select for cold/warm season
    if   season == 'all':
        d0, d1 = 0, None
    if   season == 'cold':
        d0 = (datetime(2014, 10, 1) - datetime(year0, month0, 1)).days
        d1 = (datetime(2015, 3, 31) - datetime(year0, month0, 1)).days + 1
    elif season == 'warm':
        d0 = (datetime(2015, 4, 1) - datetime(year0, month0, 1)).days
        d1 = (datetime(2015, 9, 30) - datetime(year0, month0, 1)).days + 1

    for size in sizes1:
        for en in range(ensemble):
            Pimerg[0.5, size, en] = Pimerg[0.5, size, en][d0 : d1]
            Pmrms1[0.5, size, en] = Pmrms1[0.5, size, en][d0 : d1]
    for size in sizes2:
        for en in range(ensemble):
            Ptmpa[3, size, en] = Ptmpa[3, size, en][d0 : d1]
            Pmrms2[0.5, size, en] = Pmrms2[0.5, size, en][d0 : d1]

    # calculate higher rates
    for time in times1[1:]:
        for size in sizes1:
            for en in range(ensemble):
                Pimerg[time, size, en] = accum(Pimerg[0.5, size, en], time, base = 0.5)
                Pmrms1[time, size, en] = accum(Pmrms1[0.5, size, en], time, base = 0.5)

    # calculate 3 hr rates for MRMS(TMPA)
    for size in sizes2:
        for en in range(ensemble):
            Pmrms2[3, size, en] = accum(Pmrms2[0.5, size, en], 3, base = 0.5)

    # calculate higher rates
    for time in times2[1:]:
        for size in sizes2:
            for en in range(ensemble):
                Ptmpa[time, size, en] = accum(Ptmpa[3, size, en], time, base = 3)
                Pmrms2[time, size, en] = accum(Pmrms2[3, size, en], time, base = 3)

    return Pimerg, Pmrms1, Ptmpa, Pmrms2, thresholds1, thresholds2


def compare(x, mode, minrain):
    if   mode == 'N':
        return np.ma.less(np.ma.masked_invalid(x).flatten(), minrain)
    elif mode == 'Y':
        return np.ma.greater_equal(np.ma.masked_invalid(x).flatten(), minrain)


def contingency(Pimerg, Pmrms1, Ptmpa, Pmrms2, thresholds1, thresholds2):

    contab1 = {}

    for time in times1:
        for size in sizes1:
            for en in range(ensemble):
                contab1[time, size, en] = \
                    np.array([[np.sum(compare(Pmrms1[time, size, en].flatten(), ii, 
                                              minrain = thresholds1[time, size]) * 
                                      compare(Pimerg[time, size, en].flatten(), jj, 
                                              minrain = thresholds1[time, size])) 
                               for ii in ('Y', 'N')] for jj in ('Y', 'N')])

    contab2 = {}

    for time in times2:
        for size in sizes2:
            for en in range(ensemble):
                contab2[time, size, en] = \
                    np.array([[np.sum(compare(Pmrms2[time, size, en].flatten(), ii, 
                                              minrain = thresholds2[time, size]) * 
                                      compare(Ptmpa[time, size, en].flatten(), jj, 
                                              minrain = thresholds2[time, size])) 
                               for ii in ('Y', 'N')] for jj in ('Y', 'N')])

    return contab1, contab2


def fit_mem_ols(x, y, minrain):
    '''Performs fit to the multiplicative error model using OLS. Important: x 
       is the truth while y is the variable.'''

    import statsmodels.api as sm

    th = np.ma.filled(np.ma.greater_equal(x, minrain) *
                      np.ma.greater_equal(y, minrain), False)

    if len(x[th]) == 0 or len(y[th]) == 0:
        return np.nan, np.nan, np.nan
    
    Y = np.log(y[th])
    X = sm.add_constant(np.log(x[th]))
    
    results = sm.OLS(Y, X).fit()
    
    try:
        alpha, beta = results.params
        sigma = np.std(results.resid)
    except ValueError:
        alpha, beta, sigma = np.nan, np.nan, np.nan

    return alpha, beta, sigma


def fit_aem_ols(x, y, minrain):
    '''Performs fit to the additive error model using OLS. Important: x 
       is the truth while y is the variable.'''

    import statsmodels.api as sm

    th = np.ma.filled(np.ma.greater_equal(x, minrain) *
                      np.ma.greater_equal(y, minrain), False)

    if len(x[th]) == 0 or len(y[th]) == 0:
        return np.nan, np.nan, np.nan

    Y = y[th]
    X = sm.add_constant(x[th])

    results = sm.OLS(Y, X).fit()

    try:
        a, b = results.params
        e = np.std(results.resid)
    except ValueError:
        a, b, e = np.nan, np.nan, np.nan

    return a, b, e


def fit_mem_odr(x, y, minrain):
    '''Performs fit to the multiplicative error model using ODR. Important: x 
       is the truth while y is the variable.'''

    from scipy import odr

    th = np.ma.filled(np.ma.greater_equal(x, minrain) *
                      np.ma.greater_equal(y, minrain), False)

    if len(x[th]) == 0 or len(y[th]) == 0:
        return np.nan, np.nan, np.nan

    X = np.log(x[th])
    Y = np.log(y[th])

    def f(B, x):
        return B[0] + B[1] * x

    linear = odr.Model(f)
    mydata = odr.Data(X, Y)
    myodr = odr.ODR(mydata, linear, beta0 = [1., 0])
    myoutput = myodr.run()

    return myoutput.beta[0], myoutput.beta[1], np.sqrt(myoutput.res_var)


def fit_aem_odr(x, y, minrain):
    '''Performs fit to the additive error model using OLS. Important: x 
       is the truth while y is the variable.'''

    from scipy import odr

    th = np.ma.filled(np.ma.greater_equal(x, minrain) *
                      np.ma.greater_equal(y, minrain), False)

    if len(x[th]) == 0 or len(y[th]) == 0:
        return np.nan, np.nan, np.nan

    Y = y[th]
    X = x[th]

    def f(B, x):
        return B[0] + B[1] * x

    linear = odr.Model(f)
    mydata = odr.Data(X, Y)
    myodr = odr.ODR(mydata, linear, beta0 = [1., 0])
    myoutput = myodr.run()

    return myoutput.beta[0], myoutput.beta[1], np.sqrt(myoutput.res_var)


def calc_error(x, y, error, minrain):

    '''Calculate the error of two arrays for values above the threshold.
       Important: x is the reference while y is the estimate. Choice of error 
       type ('tot', 'abs', 'rms') normalized against mean reference.'''

    x = np.array(x).flatten()
    y = np.array(y).flatten()

    th = np.ma.filled(np.ma.greater_equal(x, minrain) *
                      np.ma.greater_equal(y, minrain), False)

    if   error == 'tot':
        err = np.ma.mean(y[th] - x[th])
    elif error == 'abs':
        err = np.ma.mean(np.ma.abs(y[th] - x[th]))
    elif error == 'rms':
        err = np.ma.sqrt(np.ma.mean((y[th] - x[th]) ** 2))
    else:
        sys.exit('Error: unknown type of error.')

    return err / np.ma.mean(x[th])


def plot(f1, f2):

    for tt, time in enumerate(times1):
        plt.plot(x1, np.mean(f1[time], 0), color = cols1[tt], lw = 0.75, 
                 label = '%s h' % time)
    for tt, time in enumerate(times2):
        l, = plt.plot(x2, np.mean(f2[time], 0), color = cols2[tt], lw = 1)
        l.set_dashes(dashes)
    plt.xlabel('box length (°)')
    plt.grid()

    return None


def getMonths(year, year0, month0, year1, month1):
    if   year0 == year1: return range(month0, month1 + 1)
    elif year  == year0: return range(month0, 13)
    elif year  == year1: return range(1, month1 + 1)
    else               : return range(1, 13)


def get_CONUS_rqi():

    mrmspath = '/media/jackson/Vault/MRMS/Level3/'

    datapath = '/home/jackson/Work/Data/Scale Analysis/'
    file_rqi = '%sRQI.npy' % (datapath,)

    if os.path.exists(file_rqi):

        R = np.load(file_rqi)

    else:

        import warnings
        warnings.filterwarnings("ignore")   # ignore the empty slice warnings

        R = np.full([365, 700, 350], np.nan, dtype = np.float32)

        year = 2015
        dor = 0
        for month in range(1, 13):
            for day in range(1, monthrange(year, month)[1] + 1):

                inpath = '%s%4d/%02d/' % (mrmspath, year, month)
                ts = '%4d%02d%02d.%02d%02d00' % (year, month, day, 0, 0)
                infile = '30MRQI.%s.asc.gz' % (ts,)

                try:
                    with gzip.open('%s%s' % (inpath, infile), 'rb') as f:
                        data = f.readlines()[6:]
                    r = np.array([[float(ii) for ii in jj.split()] for jj in data])[::-1]
                    r[r < 0] = np.nan
                    R[dor] = np.array([[np.nanmean(r[ii * 10 : (ii + 1) * 10, 
                                                     jj * 10 : (jj + 1) * 10])
                                   for ii in range(350)] for jj in range(700)])
                except FileNotFoundError:
                    pass

                dor += 1

        R = np.nanmean(R, 0)
        np.save(file_rqi, R)

    return R
