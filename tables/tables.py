#!/usr/local/bin/python

# This script plots the diagrams in this paper.

import numpy as np
import sys
import pickle
from datetime import datetime
from funcs import *


#--- SORTING INPUT OPTIONS ---#


option = sys.argv[1:]

inst_options = ['inst_contab', 'inst_score']
rate_options = ['rate_corr', 'rate_error', 'rate_mem']

if   option == ['all']:
    option = inst_options + rate_options
elif option == ['inst']:
    option = inst_options
elif option == ['rate']:
    option = rate_options


#--- PRELIMINARIES ---#


# setting up the parameters
year0, year1 = 2014, 2015
month0, month1 = 4, 9
nday = (datetime(year1, month1 + 1, 1) - datetime(year0, month0, 1)).days
times1 = (0.5, 1, 3, 6, 12, 24) # IMERG times (hours) to calculate
sizes1 = list(range(1, 26))     # IMERG sizes to calculate
times2 = (3, 6, 12, 24)         # TMPA times (hours) to calculate
sizes2 = list(range(1, 11))     # TMPA sizes to calculate
x1 = np.array(sizes1) * 0.1    # convert IMERG sizes to units of degree
x2 = np.array(sizes2) * 0.25   # convert TMPA sizes to units of degree
ensemble = 100                 # ensemble size

# define the thresholds in the structure accepted by the functions
thresholds = np.hstack([1e-6, np.arange(0.05, 1.01, 0.05)])
minrain = {}
for th in thresholds:
    minrain[th] = {}
    for time in times1:
        for size in sizes1:
            minrain[th][time, size] = th


#--- RAIN INSTANCES ---#


if 'inst_contab' in option:

    if 'Pimerg' not in locals(): 
        Pimerg, Pmrms1, Ptmpa, Pmrms2 = read_data()

    fhit = lambda x : x[0, 0] / np.sum(x) * 100
    fmiss = lambda x : x[1, 0] / np.sum(x) * 100
    ffalse = lambda x : x[0, 1] / np.sum(x) * 100
    fcorrej = lambda x : x[1, 1] / np.sum(x) * 100

    for th in thresholds:

        contab1, contab2 = contingency(Pimerg, Pmrms1, Ptmpa, Pmrms2, 
                                       minrain[th], minrain[th])

        hit1, miss1, false1, correj1 = {}, {}, {}, {}

        for time in times1:
            hit1[time] = np.array([[fhit(contab1[time, size, en])
                                  for size in sizes1] for en in range(ensemble)])
            miss1[time] = np.array([[fmiss(contab1[time, size, en])
                                   for size in sizes1] for en in range(ensemble)])
            false1[time] = np.array([[ffalse(contab1[time, size, en])
                                    for size in sizes1] for en in range(ensemble)])
            correj1[time] = np.array([[fcorrej(contab1[time, size, en])
                                     for size in sizes1] for en in range(ensemble)])

        hit2, miss2, false2, correj2 = {}, {}, {}, {}

        for time in times2:
            hit2[time] = np.array([[fhit(contab2[time, size, en])
                                  for size in sizes2] for en in range(ensemble)])
            miss2[time] = np.array([[fmiss(contab2[time, size, en])
                                   for size in sizes2] for en in range(ensemble)])
            false2[time] = np.array([[ffalse(contab2[time, size, en])
                                    for size in sizes2] for en in range(ensemble)])
            correj2[time] = np.array([[fcorrej(contab2[time, size, en])
                                     for size in sizes2] for en in range(ensemble)])

        np.savetxt('output/imerg_hit_%4.2f.txt' % th, 
                   np.array([np.mean(hit1[time], 0) for time in times1]).T,
                   fmt = '%8.3f')
        np.savetxt('output/imerg_miss_%4.2f.txt' % th, 
                   np.array([np.mean(miss1[time], 0) for time in times1]).T,
                   fmt = '%8.3f')
        np.savetxt('output/imerg_false_%4.2f.txt' % th, 
                   np.array([np.mean(false1[time], 0) for time in times1]).T,
                   fmt = '%8.3f')
        np.savetxt('output/imerg_corneg_%4.2f.txt' % th, 
                   np.array([np.mean(correj1[time], 0) for time in times1]).T,
                   fmt = '%8.3f')
        np.savetxt('output/tmpa_hit_%4.2f.txt' % th, 
                   np.array([np.mean(hit2[time], 0) for time in times2]).T,
                   fmt = '%8.3f')
        np.savetxt('output/tmpa_miss_%4.2f.txt' % th, 
                   np.array([np.mean(miss2[time], 0) for time in times2]).T,
                   fmt = '%8.3f')
        np.savetxt('output/tmpa_false_%4.2f.txt' % th, 
                   np.array([np.mean(false2[time], 0) for time in times2]).T,
                   fmt = '%8.3f')
        np.savetxt('output/tmpa_corneg_%4.2f.txt' % th, 
                   np.array([np.mean(correj2[time], 0) for time in times2]).T,
                   fmt = '%8.3f')


if 'inst_score' in option:

    if 'Pimerg' not in locals(): 
        Pimerg, Pmrms1, Ptmpa, Pmrms2 = read_data()

    # probability of detection
    pod = lambda x : (x[0, 0] / (x[0, 0] + x[1, 0]))
    far = lambda x : (x[0, 1] / (x[0, 0] + x[0, 1]))
    bid = lambda x : ((x[0, 0] + x[0, 1]) / (x[0, 0] + x[1, 0]))
    def hss(x):
        N = np.sum(x)
        exp = 1 / N * ((x[0, 0] + x[1, 0]) * (x[0, 0] + x[0, 1]) + 
                       (x[1, 1] + x[1, 0]) * (x[1, 1] + x[0, 1]))
        return ((x[0, 0] + x[1, 1] - exp) / (N - exp))

    for th in thresholds:

        contab1, contab2 = contingency(Pimerg, Pmrms1, Ptmpa, Pmrms2, 
                                       minrain[th], minrain[th])

        pod1 = {}
        for time in times1:
            pod1[time] = np.array([[pod(contab1[time, size, en])
                                    for size in sizes1] for en in range(ensemble)])

        pod2 = {}
        for time in times2:
            pod2[time] = np.array([[pod(contab2[time, size, en])
                                    for size in sizes2] for en in range(ensemble)])

        far1 = {}
        for time in times1:
            far1[time] = np.array([[far(contab1[time, size, en])
                                    for size in sizes1] for en in range(ensemble)])

        far2 = {}
        for time in times2:
            far2[time] = np.array([[far(contab2[time, size, en])
                                   for size in sizes2] for en in range(ensemble)])

        bid1 = {}
        for time in times1:
            bid1[time] = np.array([[bid(contab1[time, size, en])
                                    for size in sizes1] for en in range(ensemble)])

        bid2 = {}
        for time in times2:
            bid2[time] = np.array([[bid(contab2[time, size, en])
                                    for size in sizes2] for en in range(ensemble)])

        hss1 = {}
        for time in times1:
            hss1[time] = np.array([[hss(contab1[time, size, en])
                                    for size in sizes1] for en in range(ensemble)])

        hss2 = {}
        for time in times2:
            hss2[time] = np.array([[hss(contab2[time, size, en])
                                    for size in sizes2] for en in range(ensemble)])

        np.savetxt('output/imerg_pod_%4.2f.txt' % th, 
                   np.array([np.mean(pod1[time], 0) for time in times1]).T,
                   fmt = '%8.3f')
        np.savetxt('output/imerg_far_%4.2f.txt' % th, 
                   np.array([np.mean(far1[time], 0) for time in times1]).T,
                   fmt = '%8.3f')
        np.savetxt('output/imerg_bid_%4.2f.txt' % th, 
                   np.array([np.mean(bid1[time], 0) for time in times1]).T,
                   fmt = '%8.3f')
        np.savetxt('output/imerg_hss_%4.2f.txt' % th, 
                   np.array([np.mean(hss1[time], 0) for time in times1]).T,
                   fmt = '%8.3f')
        np.savetxt('output/tmpa_pod_%4.2f.txt' % th, 
                   np.array([np.mean(pod2[time], 0) for time in times2]).T,
                   fmt = '%8.3f')
        np.savetxt('output/tmpa_far_%4.2f.txt' % th, 
                   np.array([np.mean(far2[time], 0) for time in times2]).T,
                   fmt = '%8.3f')
        np.savetxt('output/tmpa_bid_%4.2f.txt' % th, 
                   np.array([np.mean(bid2[time], 0) for time in times2]).T,
                   fmt = '%8.3f')
        np.savetxt('output/tmpa_hss_%4.2f.txt' % th, 
                   np.array([np.mean(hss2[time], 0) for time in times2]).T,
                   fmt = '%8.3f')



#--- DIAGRAMS ON RAIN RATES ---#


if 'rate_corr' in option:

    if 'Pimerg' not in locals(): 
        Pimerg, Pmrms1, Ptmpa, Pmrms2 = read_data()

    def corrcoef(x, y, minrain):
        x = x.flatten()
        y = y.flatten()
        th = np.ma.filled((x >= minrain) * (y >= minrain), False)
        return np.ma.corrcoef(x[th], y[th])[0, 1]

    for th in thresholds:

        corrs1, corrs2 = {}, {}

        for time in times1:
            corrs1[time] = np.array([[corrcoef(Pmrms1[time, size, en], Pimerg[time, size, en], 
                                      minrain[th][time, size]) for size in sizes1] 
                                     for en in range(ensemble)])

        for time in times2:
            corrs2[time] = np.array([[corrcoef(Pmrms2[time, size, en], Ptmpa[time, size, en], 
                                      minrain[th][time, size]) for size in sizes2] 
                                     for en in range(ensemble)])

        np.savetxt('output/imerg_corr_%4.2f.txt' % th, 
                   np.array([np.mean(corrs1[time], 0) for time in times1]).T,
                   fmt = '%8.3f')
        np.savetxt('output/tmpa_corr_%4.2f.txt' % th, 
                   np.array([np.mean(corrs2[time], 0) for time in times2]).T,
                   fmt = '%8.3f')


if 'rate_error' in option:

    if 'Pimerg' not in locals(): 
        Pimerg, Pmrms1, Ptmpa, Pmrms2 = read_data()

    for th in thresholds:

        errtot1, errtot2 = {}, {}

        for time in times1:
            errtot1[time] = np.array([[calc_error(Pmrms1[time, size, en].flatten(),
                                                  Pimerg[time, size, en].flatten(),
                                                  error = 'tot',
                                                  minrain = minrain[th][time, size])
                                     for size in sizes1] for en in range(ensemble)])

        for time in times2:
            errtot2[time] = np.array([[calc_error(Pmrms2[time, size, en].flatten(), 
                                                  Ptmpa[time, size, en].flatten(),
                                                  error = 'tot',
                                                  minrain = minrain[th][time, size])
                                     for size in sizes2] for en in range(ensemble)])

        errabs1, errabs2 = {}, {}

        for time in times1:
            errabs1[time] = np.array([[calc_error(Pmrms1[time, size, en].flatten(),
                                                  Pimerg[time, size, en].flatten(),
                                                  error = 'abs',
                                                  minrain = minrain[th][time, size])
                                     for size in sizes1] for en in range(ensemble)])

        for time in times2:
            errabs2[time] = np.array([[calc_error(Pmrms2[time, size, en].flatten(), 
                                                  Ptmpa[time, size, en].flatten(),
                                                  error = 'abs',
                                                  minrain = minrain[th][time, size])
                                     for size in sizes2] for en in range(ensemble)])

        errrms1, errrms2 = {}, {}

        for time in times1:
            errrms1[time] = np.array([[calc_error(Pmrms1[time, size, en].flatten(),
                                                  Pimerg[time, size, en].flatten(),
                                                  error = 'rms',
                                                  minrain = minrain[th][time, size])
                                     for size in sizes1] for en in range(ensemble)])

        for time in times2:
            errrms2[time] = np.array([[calc_error(Pmrms2[time, size, en].flatten(), 
                                                  Ptmpa[time, size, en].flatten(),
                                                  error = 'rms',
                                                  minrain = minrain[th][time, size])
                                     for size in sizes2] for en in range(ensemble)])

        np.savetxt('output/imerg_nme_%4.2f.txt' % th, 
                   np.array([np.mean(errtot1[time], 0) for time in times1]).T,
                   fmt = '%8.3f')
        np.savetxt('output/imerg_nmae_%4.2f.txt' % th, 
                   np.array([np.mean(errabs1[time], 0) for time in times1]).T,
                   fmt = '%8.3f')
        np.savetxt('output/imerg_nrmse_%4.2f.txt' % th, 
                   np.array([np.mean(errrms1[time], 0) for time in times1]).T,
                   fmt = '%8.3f')
        np.savetxt('output/tmpa_nme_%4.2f.txt' % th, 
                   np.array([np.mean(errtot2[time], 0) for time in times2]).T,
                   fmt = '%8.3f')
        np.savetxt('output/tmpa_nmae_%4.2f.txt' % th, 
                   np.array([np.mean(errabs2[time], 0) for time in times2]).T,
                   fmt = '%8.3f')
        np.savetxt('output/tmpa_nrmse_%4.2f.txt' % th, 
                   np.array([np.mean(errrms2[time], 0) for time in times2]).T,
                   fmt = '%8.3f')


if 'rate_mem' in option:

    if 'Pimerg' not in locals(): 
        Pimerg, Pmrms1, Ptmpa, Pmrms2 = read_data()

    for th in thresholds:

        alpha1, beta1, sigma1 = {}, {}, {}

        fit_mem = fit_mem_ols    # choose OLS instead of ODR

        for time in times1:
            params = [[fit_mem(Pmrms1[time, size, en], Pimerg[time, size, en], 
                               minrain = minrain[th][time, size])
                        for size in sizes1] for en in range(ensemble)]
            alpha1[time] = np.array([[jj[0] for jj in ii] for ii in params])
            beta1[time] = np.array([[jj[1] for jj in ii] for ii in params])
            sigma1[time] = np.array([[jj[2] for jj in ii] for ii in params])

        alpha2, beta2, sigma2 = {}, {}, {}

        for time in times2:
            params = [[fit_mem(Pmrms2[time, size, en], Ptmpa[time, size, en],
                               minrain = minrain[th][time, size])
                        for size in sizes2] for en in range(ensemble)]
            alpha2[time] = np.array([[jj[0] for jj in ii] for ii in params])
            beta2[time] = np.array([[jj[1] for jj in ii] for ii in params])
            sigma2[time] = np.array([[jj[2] for jj in ii] for ii in params])

        np.savetxt('output/imerg_mem-alpha_%4.2f.txt' % th, 
                   np.array([np.mean(alpha1[time], 0) for time in times1]).T,
                   fmt = '%8.3f')
        np.savetxt('output/imerg_mem-beta_%4.2f.txt' % th, 
                   np.array([np.mean(beta1[time], 0) for time in times1]).T,
                   fmt = '%8.3f')
        np.savetxt('output/imerg_mem-sigma_%4.2f.txt' % th, 
                   np.array([np.mean(sigma1[time], 0) for time in times1]).T,
                   fmt = '%8.3f')
        np.savetxt('output/tmpa_mem-alpha_%4.2f.txt' % th, 
                   np.array([np.mean(alpha2[time], 0) for time in times2]).T,
                   fmt = '%8.3f')
        np.savetxt('output/tmpa_mem-beta_%4.2f.txt' % th, 
                   np.array([np.mean(beta2[time], 0) for time in times2]).T,
                   fmt = '%8.3f')
        np.savetxt('output/tmpa_mem-sigma_%4.2f.txt' % th, 
                   np.array([np.mean(sigma2[time], 0) for time in times2]).T,
                   fmt = '%8.3f')
