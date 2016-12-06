#!/usr/local/bin/python

# This script plots the diagrams in this paper.

import numpy as np
import matplotlib.pyplot as plt
import brewer2mpl
import sys
import pickle
from calendar import monthrange
from datetime import datetime
from funcs import *


#--- SORTING INPUT OPTIONS ---#


option = sys.argv[1:]

bkgd_options = ['bkgd_map', 'bkgd_example', 'bkgd_fraction', 'bkgd_thresholds']
inst_options = ['inst_contab', 'inst_score']
rate_options = ['rate_corr', 'rate_error', 'rate_mem']
sup_options  = ['sup_uncal', 'sup_mem']

if   option == ['all']:
    option = bkgd_options + inst_options + rate_options + sup_options
elif option == ['bkgd']:
    option = bkgd_options
elif option == ['inst']:
    option = inst_options
elif option == ['rate']:
    option = rate_options
elif option == ['sup']:
    option = sup_options


#--- PRELIMINARIES ---#


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


#--- DIAGRAMS ON BACKGROUND ---#


if 'bkgd_map' in option:

    from mpl_toolkits.basemap import Basemap
    from matplotlib.patches import Rectangle

    m = Basemap(projection = 'cyl', llcrnrlat = 20, urcrnrlat = 60, 
                llcrnrlon = -130, urcrnrlon = -70, lat_ts = 40, 
                resolution = 'i')

    plt.figure(figsize = (scol, 0.75 * scol))

    lonedges = np.arange(-130, -59.99, 0.1)
    latedges = np.arange(  20,  55.01, 0.1)
    R = get_CONUS_rqi()
    plt.pcolormesh(lonedges, latedges, np.ma.masked_invalid(R).T, 
                   vmin = 0, vmax = 100, cmap = plt.cm.viridis_r, 
                   rasterized = True)
    m.drawcountries(linewidth = 0.5, color = 'k', zorder = 1)
    m.drawstates(linewidth = 0.5, color = '0.5', zorder = 1)
    m.drawcoastlines(linewidth = 0.5, zorder = 1)
    plt.gca().add_patch(Rectangle((-93.5, 30), 10, 11.5, facecolor = 'none', 
                        linewidth = 1, edgecolor = 'r', zorder = 2))
    plt.text(-94.5, 30, '30°N', size = 8, color = 'r', 
             ha = 'right', va = 'center')
    plt.text(-94.5, 41.5, '41.5°N', size = 8, color = 'r', 
             ha = 'right', va = 'center')
    plt.text(-93.5, 29, '93.5°W', size = 8, color = 'r', 
             ha = 'center', va = 'top')
    plt.text(-83.5, 29, '83.5°W', size = 8, color = 'r', 
             ha = 'center', va = 'top')
    m.drawparallels(np.arange(20, 61, 20), linewidth = 0, 
                    labels = [1, 0, 0, 0])
    m.drawmeridians(np.arange(-130, -69, 20), linewidth = 0, 
                    labels = [0, 0, 0, 1])

    plt.colorbar(cax = plt.axes([0.95, 0.15, 0.025, 0.7]))

    plt.savefig('fig_bkgd_map.pdf')
    plt.close()


if 'bkgd_example' in option:

    if 'Pimerg' not in locals(): 
        Pimerg, Pmrms1, Ptmpa, Pmrms2, thresholds1, thresholds2 = read_data()

    example_x1 = np.ma.masked_array([ii for en in range(10) 
                             for ii in Pimerg[0.5, 1, en].flatten()])
    example_y1 = np.ma.masked_array([ii for en in range(10) 
                             for ii in Pmrms1[0.5, 1, en].flatten()])
    example_x2 = np.ma.masked_array([ii for en in range(60) 
                             for ii in Pimerg[24, 1, en].flatten()])
    example_y2 = np.ma.masked_array([ii for en in range(60) 
                             for ii in Pmrms1[24, 1, en].flatten()])
    example_x3 = np.ma.masked_array([ii for en in range(10) 
                             for ii in Pimerg[0.5, 25, en].flatten()])
    example_y3 = np.ma.masked_array([ii for en in range(10) 
                             for ii in Pmrms1[0.5, 25, en].flatten()])
    example_x4 = np.ma.masked_array([ii for en in range(60) 
                             for ii in Pimerg[24, 25, en].flatten()])
    example_y4 = np.ma.masked_array([ii for en in range(60) 
                             for ii in Pmrms1[24, 25, en].flatten()])

    def scatter_log(x, y, size = 1, threshold = 0.1):

        th = np.ma.filled(np.ma.greater_equal(x, threshold) *
                          np.ma.greater_equal(y, threshold), False)

        plt.scatter(x[th], y[th], s = size, color = 'k', edgecolor = 'none',
                    rasterized = True)
        plt.plot([0.1, 100], [0.1, 100], color = '0.5', ls = '--')
        plt.gca().set_yscale('log')
        plt.gca().set_xscale('log')
        plt.axis([0.1, 100, 0.1, 100])
        plt.grid()

        return None

    plt.figure(figsize = (scol, scol))
    plt.subplots_adjust(hspace = 0.3, wspace = 0.3)
    plt.subplot(221)
    scatter_log(example_x1, example_y1)
    plt.text(0.02, 0.96, '(a)', ha = 'left', va = 'top',
             transform = plt.gca().transAxes)
    plt.ylabel('IMERG (mm / h)')
    plt.title('0.1°, 0.5 h')
    plt.subplot(222)
    scatter_log(example_x2, example_y2)
    plt.text(0.02, 0.96, '(b)', ha = 'left', va = 'top',
             transform = plt.gca().transAxes)
    plt.title('0.1°, 24 h')
    plt.subplot(223)
    scatter_log(example_x3, example_y3)
    plt.text(0.02, 0.96, '(c)', ha = 'left', va = 'top',
             transform = plt.gca().transAxes)
    plt.xlabel('reference (mm / h)')
    plt.ylabel('IMERG (mm / h)')
    plt.title('2.5°, 0.5 h')
    plt.subplot(224)
    scatter_log(example_x4, example_y4)
    plt.text(0.02, 0.96, '(d)', ha = 'left', va = 'top',
             transform = plt.gca().transAxes)
    plt.xlabel('reference (mm / h)')
    plt.title('2.5°, 24 h')
    plt.savefig('fig_bkgd_example.pdf')
    plt.close()


if 'bkgd_fraction' in option:

    if 'Pimerg' not in locals(): 
        Pimerg, Pmrms1, Ptmpa, Pmrms2, thresholds1, thresholds2 = read_data()

    def count_rain_frac(X, time, minrain = 0.2):
        return np.sum(compare(X, 'Y', minrain = minrain)) / (nday * 24 / time)

    frac1 = {}
    for time in times1:
        for size in sizes1:
            frac1[time, size] = np.mean([count_rain_frac(Pmrms1[time, size, en], 
                                         time, minrain = 0.2) 
                                         for en in range(ensemble)])

    plt.figure(figsize = (scol, 0.75 * scol))
    for tt, time in enumerate(times1):
        plt.plot(x1, [frac1[time, size] for size in sizes1], color = cols1[tt], 
                 lw = 0.75, label = '%s h' % time)
    plt.xlabel('box length (°)')
    plt.ylabel('fraction of events with at least 0.2 mm / h')
    plt.grid()
    plt.legend(loc = 4, ncol = 2)
    plt.savefig('fig_bkgd_fraction.pdf')
    plt.close()


if 'bkgd_thresholds' in option:

    minrain1 = 0.2         # threshold for the highest resolution of IMERG
    minrain2 = 0.2         # threshold for the highest resolution of TMPA

    thresholds1 = {}
    for time in times1:
        for size in sizes1:
            thresholds1[time, size] = minrain1 / (size * np.sqrt(time / 0.5))

    thresholds2 = {}
    for time in times2:
        for size in sizes2:
            thresholds2[time, size] = minrain2 / (size * 2.5 * np.sqrt(time / 0.5))

    plt.figure(figsize = (scol, 0.75 * scol))
    for tt, time in enumerate(times1):
        plt.plot(x1, [thresholds1[time, size] for size in sizes1], 
                 color = cols1[tt], lw = 0.75, label = '%s h' % time)
    for tt, time in enumerate(times2):
        l, = plt.plot(x2, [thresholds2[time, size] for size in sizes2],
                 color = cols2[tt], lw = 1)
        l.set_dashes(dashes)
    plt.xlabel('box length (°)')
    plt.ylabel('threshold (mm / h)')
    plt.legend(loc = 1, ncol = 2)
    plt.grid()
    plt.savefig('fig_bkgd_thresholds.pdf')
    plt.close()

    pickle.dump(thresholds1, open('data_thresholds1.p', 'wb'))
    pickle.dump(thresholds2, open('data_thresholds2.p', 'wb'))


if 'bkgd_thresholds_hss' in option:

    from signal_smooth import smooth

    if 'Pimerg' not in locals(): 
        Pimerg, Pmrms1, Ptmpa, Pmrms2, _, _ = read_data()

    def hss(x):
        x = np.ma.masked_invalid(x)
        N = np.sum(x)
        exp = 1 / N * ((x[0, 0] + x[1, 0]) * (x[0, 0] + x[0, 1]) + 
                       (x[1, 1] + x[1, 0]) * (x[1, 1] + x[0, 1]))
        return ((x[0, 0] + x[1, 1] - exp) / (N - exp))

    minrains = np.arange(0.01, 0.5001, 0.01)

    thresholds1 = {}
    for time in times1:
        for size in sizes1:

            hss1 = {}
            for minrain in minrains:
                C = np.array([[[np.sum(compare(Pmrms1[time, size, en],
                                               ii, minrain = minrain) * 
                                       compare(Pimerg[time, size, en], 
                                               jj, minrain = minrain)) 
                                for ii in ('Y', 'N')] for jj in ('Y', 'N')]
                              for en in range(ensemble)])
                hss1[minrain] = hss(np.sum(C, 0))

            thresholds1[time, size] = max(hss1, key = hss1.get)


    thresholds2 = {}
    for time in times2:
        for size in sizes2:

            hss2 = {}
            for minrain in minrains:
                C = np.array([[[np.sum(compare(Pmrms2[time, size, en],
                                               ii, minrain = minrain) * 
                                       compare(Ptmpa[time, size, en], 
                                               jj, minrain = minrain)) 
                                for ii in ('Y', 'N')] for jj in ('Y', 'N')]
                              for en in range(ensemble)])
                hss2[minrain] = hss(np.sum(C, 0))

            thresholds2[time, size] = max(hss2, key = hss2.get)

    thresholds_smoothed1 = {}
    for time in times1:
        y1 = smooth(np.array([thresholds1[time, size] for size in sizes1]), 11, window = 'flat')
        for ss, size in enumerate(sizes1):
            thresholds_smoothed1[time, size] = y1[ss]
    thresholds_smoothed2 = {}
    for time in times2:
        y2 = smooth(np.array([thresholds1[time, size] for size in sizes2]), 3, window = 'flat')
        for ss, size in enumerate(sizes2):
            thresholds_smoothed2[time, size] = y2[ss]

    plt.figure(figsize = (scol, 0.75 * scol))
    for tt, time in enumerate(times1):
        plt.plot(x1, [thresholds_smoothed1[time, size] for size in sizes1],
                 color = cols1[tt], lw = 0.75, label = '%s h' % time)
    for tt, time in enumerate(times2):
        l, = plt.plot(x2, [thresholds_smoothed2[time, size] for size in sizes2],
                 color = cols2[tt], lw = 1)
        l.set_dashes(dashes)
    plt.xlabel('box length (°)')
    plt.ylabel('threshold (mm / h)')
    plt.legend(loc = 2, ncol = 2, framealpha = 0.5)
    plt.grid()
    plt.savefig('fig_bkgd_thresholds_hss.pdf')
    plt.close()

    pickle.dump(thresholds1, open('data_thresholds1_hss.p', 'wb'))
    pickle.dump(thresholds2, open('data_thresholds2_hss.p', 'wb'))


#--- DIAGRAMS ON RAIN INSTANCES ---#


if 'inst_contab' in option:

    if 'Pimerg' not in locals(): 
        Pimerg, Pmrms1, Ptmpa, Pmrms2, thresholds1, thresholds2 = read_data()
    if 'contab1' not in locals():
        contab1, contab2 = contingency(Pimerg, Pmrms1, Ptmpa, Pmrms2, 
                                       thresholds1, thresholds2)

    fhit = lambda x : x[0, 0] / np.sum(x) * 100
    fmiss = lambda x : x[1, 0] / np.sum(x) * 100
    ffalse = lambda x : x[0, 1] / np.sum(x) * 100
    fcorrej = lambda x : x[1, 1] / np.sum(x) * 100

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

    plt.figure(figsize = (dcol, 0.75 * dcol))
    plt.subplot(221)
    plot(hit1, hit2)
    plt.setp(plt.gca(), 'xticklabels', [])
    plt.xlabel('')
    plt.ylabel('hit (%)')
    plt.subplot(222)
    plot(miss1, miss2)
    plt.setp(plt.gca(), 'xticklabels', [])
    plt.xlabel('')
    plt.ylabel('miss (%)')
    plt.legend(loc = 4, ncol = 2)
    plt.subplot(223)
    plot(false1, false2)
    plt.ylabel('false alarm (%)')
    plt.subplot(224)
    plot(correj1, correj2)
    plt.ylabel('correct rejection (%)')
    plt.savefig('fig_inst_contab.pdf')
    plt.close()


if 'inst_score' in option:

    if 'Pimerg' not in locals(): 
        Pimerg, Pmrms1, Ptmpa, Pmrms2, thresholds1, thresholds2 = read_data()
    if 'contab1' not in locals():
        contab1, contab2 = contingency(Pimerg, Pmrms1, Ptmpa, Pmrms2, 
                                       thresholds1, thresholds2)

    # probability of detection
    pod = lambda x : (x[0, 0] / (x[0, 0] + x[1, 0]))
    far = lambda x : (x[0, 1] / (x[0, 0] + x[0, 1]))
    bid = lambda x : ((x[0, 0] + x[0, 1]) / (x[0, 0] + x[1, 0]))
    def hss(x):
        N = np.sum(x)
        exp = 1 / N * ((x[0, 0] + x[1, 0]) * (x[0, 0] + x[0, 1]) + 
                       (x[1, 1] + x[1, 0]) * (x[1, 1] + x[0, 1]))
        return ((x[0, 0] + x[1, 1] - exp) / (N - exp))

    pod1, far1, bid1, hss1 = {}, {}, {}, {}
    for time in times1:
        pod1[time] = np.array([[pod(contab1[time, size, en])
                                for size in sizes1] for en in range(ensemble)])
        far1[time] = np.array([[far(contab1[time, size, en])
                                for size in sizes1] for en in range(ensemble)])
        bid1[time] = np.array([[bid(contab1[time, size, en])
                                for size in sizes1] for en in range(ensemble)])
        hss1[time] = np.array([[hss(contab1[time, size, en])
                                for size in sizes1] for en in range(ensemble)])

    pod2, far2, bid2, hss2 = {}, {}, {}, {}
    for time in times2:
        pod2[time] = np.array([[pod(contab2[time, size, en])
                                for size in sizes2] for en in range(ensemble)])
        far2[time] = np.array([[far(contab2[time, size, en])
                               for size in sizes2] for en in range(ensemble)])
        bid2[time] = np.array([[bid(contab2[time, size, en])
                                for size in sizes2] for en in range(ensemble)])
        hss2[time] = np.array([[hss(contab2[time, size, en])
                                for size in sizes2] for en in range(ensemble)])

    plt.figure(figsize = (dcol, 0.75 * dcol))
    plt.subplots_adjust(wspace = 0.25)
    plt.subplot(221)
    plot(pod1, pod2)
    plt.setp(plt.gca(), 'xticklabels', [])
    plt.xlabel('')
    plt.ylabel('probability of detection')
    plt.subplot(222)
    plot(far1, far2)
    plt.setp(plt.gca(), 'xticklabels', [])
    plt.xlabel('')
    plt.ylabel('false alarm ratio')
    plt.legend(loc = 1, ncol = 2)
    plt.subplot(223)
    plot(bid1, bid2)
    plt.ylabel('bias in detection')
    plt.subplot(224)
    plot(hss1, hss2)
    plt.ylabel('Heidke skill score')
    plt.savefig('fig_inst_score.pdf')
    plt.close()


#--- DIAGRAMS ON RAIN RATES ---#


if 'rate_corr' in option:

    if 'Pimerg' not in locals(): 
        Pimerg, Pmrms1, Ptmpa, Pmrms2, thresholds1, thresholds2 = read_data()

    def corrcoef(x, y, minrain):
        x = np.nan_to_num(x).flatten()
        y = np.nan_to_num(y).flatten()
        th = np.ma.filled((x >= minrain) * (y >= minrain), False)
        return np.ma.corrcoef(x[th], y[th])[0, 1]

    corrs1, corrs2 = {}, {}

    for time in times1:
        corrs1[time] = np.array([[corrcoef(Pmrms1[time, size, en], Pimerg[time, size, en], 
                                  thresholds1[time, size]) for size in sizes1] 
                                 for en in range(ensemble)])

    for time in times2:
        corrs2[time] = np.array([[corrcoef(Pmrms2[time, size, en], Ptmpa[time, size, en], 
                                  thresholds2[time, size]) for size in sizes2] 
                                 for en in range(ensemble)])

    plt.figure(figsize = (scol, 0.75 * scol))
    plot(corrs1, corrs2)
    plt.ylabel('correlation')
    plt.legend(loc = 4)
    plt.savefig('fig_rate_corr.pdf')
    plt.close()


if 'rate_error' in option:

    if 'Pimerg' not in locals(): 
        Pimerg, Pmrms1, Ptmpa, Pmrms2, thresholds1, thresholds2 = read_data()

    errtot1, errabs1, errrms1 = {}, {}, {}

    for time in times1:
        errtot1[time] = np.array([[calc_error(Pmrms1[time, size, en].flatten(),
                                              Pimerg[time, size, en].flatten(),
                                              error = 'tot',
                                              minrain = thresholds1[time, size])
                                 for size in sizes1] for en in range(ensemble)])
        errabs1[time] = np.array([[calc_error(Pmrms1[time, size, en].flatten(),
                                              Pimerg[time, size, en].flatten(),
                                              error = 'abs',
                                              minrain = thresholds1[time, size])
                                 for size in sizes1] for en in range(ensemble)])
        errrms1[time] = np.array([[calc_error(Pmrms1[time, size, en].flatten(),
                                              Pimerg[time, size, en].flatten(),
                                              error = 'rms',
                                              minrain = thresholds1[time, size])
                                 for size in sizes1] for en in range(ensemble)])

    errtot2, errabs2, errrms2 = {}, {}, {}

    for time in times2:
        errtot2[time] = np.array([[calc_error(Pmrms2[time, size, en].flatten(), 
                                              Ptmpa[time, size, en].flatten(),
                                              error = 'tot',
                                              minrain = thresholds2[time, size])
                                 for size in sizes2] for en in range(ensemble)])
        errabs2[time] = np.array([[calc_error(Pmrms2[time, size, en].flatten(), 
                                              Ptmpa[time, size, en].flatten(),
                                              error = 'abs',
                                              minrain = thresholds1[time, size])
                                 for size in sizes2] for en in range(ensemble)])
        errrms2[time] = np.array([[calc_error(Pmrms2[time, size, en].flatten(), 
                                              Ptmpa[time, size, en].flatten(),
                                              error = 'rms',
                                              minrain = thresholds1[time, size])
                                 for size in sizes2] for en in range(ensemble)])

    plt.figure(figsize = (scol, 1.33 * scol))
    plt.subplot(311)
    plot(errtot1, errtot2)
    plt.setp(plt.gca(), 'xticklabels', [])
    plt.xlabel('')
    plt.ylabel('norm. mean error')
    plt.subplot(312)
    plot(errabs1, errabs2)
    plt.setp(plt.gca(), 'xticklabels', [])
    plt.xlabel('')
    plt.ylabel('norm. mean abs. error')
    plt.legend(loc = 1, ncol = 2)
    plt.subplot(313)
    plot(errrms1, errrms2)
    plt.ylabel('norm. RMSE')
    plt.savefig('fig_rate_error.pdf')
    plt.close()


if 'rate_mem' in option:

    if 'Pimerg' not in locals(): 
        Pimerg, Pmrms1, Ptmpa, Pmrms2, thresholds1, thresholds2 = read_data()

    alpha1, beta1, sigma1 = {}, {}, {}

    fit_mem = fit_mem_ols    # choose OLS instead of ODR

    for time in times1:
        params = [[fit_mem(Pmrms1[time, size, en], Pimerg[time, size, en], 
                           minrain = thresholds1[time, size])
                    for size in sizes1] for en in range(ensemble)]
        alpha1[time] = np.array([[jj[0] for jj in ii] for ii in params])
        beta1[time] = np.array([[jj[1] for jj in ii] for ii in params])
        sigma1[time] = np.array([[jj[2] for jj in ii] for ii in params])

    alpha2, beta2, sigma2 = {}, {}, {}

    for time in times2:
        params = [[fit_mem(Pmrms2[time, size, en], Ptmpa[time, size, en],
                           minrain = thresholds2[time, size])
                    for size in sizes2] for en in range(ensemble)]
        alpha2[time] = np.array([[jj[0] for jj in ii] for ii in params])
        beta2[time] = np.array([[jj[1] for jj in ii] for ii in params])
        sigma2[time] = np.array([[jj[2] for jj in ii] for ii in params])

    plt.figure(figsize = (scol, 1.33 * scol))

    plt.subplot(311)
    plot(alpha1, alpha2)
    plt.setp(plt.gca(), 'xticklabels', [])
    plt.legend(loc = 1, ncol = 2)
    plt.xlabel('')
    plt.ylabel(r'$\alpha$')

    plt.subplot(312)
    plot(beta1, beta2)
    plt.setp(plt.gca(), 'xticklabels', [])
    plt.xlabel('')
    plt.ylabel(r'$\beta$')

    plt.subplot(313)
    plot(sigma1, sigma2)
    plt.ylabel(r'$\sigma$')

    plt.savefig('fig_rate_mem.pdf')
    plt.close()


if 'rate_aem' in option:

    if 'Pimerg' not in locals(): 
        Pimerg, Pmrms1, Ptmpa, Pmrms2, thresholds1, thresholds2 = read_data()

    alpha1, beta1, sigma1 = {}, {}, {}

    fit_mem = fit_aem_ols

    for time in times1:
        params = [[fit_mem(Pmrms1[time, size, en], Pimerg[time, size, en], 
                           minrain = thresholds1[time, size])
                    for size in sizes1] for en in range(ensemble)]
        alpha1[time] = np.array([[jj[0] for jj in ii] for ii in params])
        beta1[time] = np.array([[jj[1] for jj in ii] for ii in params])
        sigma1[time] = np.array([[jj[2] for jj in ii] for ii in params])

    alpha2, beta2, sigma2 = {}, {}, {}

    for time in times2:
        params = [[fit_mem(Pmrms2[time, size, en], Ptmpa[time, size, en],
                           minrain = thresholds1[time, size])
                    for size in sizes2] for en in range(ensemble)]
        alpha2[time] = np.array([[jj[0] for jj in ii] for ii in params])
        beta2[time] = np.array([[jj[1] for jj in ii] for ii in params])
        sigma2[time] = np.array([[jj[2] for jj in ii] for ii in params])

    plt.figure(figsize = (scol, 1.33 * scol))

    plt.subplot(311)
    plot(alpha1, alpha2)
    plt.setp(plt.gca(), 'xticklabels', [])
    plt.legend(loc = 1, ncol = 2)
    plt.xlabel('')
    plt.ylabel('$a$')

    plt.subplot(312)
    plot(beta1, beta2)
    plt.setp(plt.gca(), 'xticklabels', [])
    plt.xlabel('')
    plt.ylabel('$b$')

    plt.subplot(313)
    plot(sigma1, sigma2)
    plt.ylabel(r'$\sigma$')

    plt.savefig('fig_rate_aem.pdf')
    plt.close()


#--- APPENDIX DIAGRAMS ---#


if 'app_mem' in option:

    memcols = ('#e41a1c', '#ff7f00', 'k', '#377eb8', '#984ea3')

    x = np.linspace(0, 100, 101)
    alphas = np.arange(-1, 1.01, 0.5)
    betas = np.arange(0.5, 1.51, 0.25)

    plt.figure(figsize = (dcol, 0.8 * scol))
    plt.subplots_adjust(wspace = 0.05)

    plt.subplot(121)
    for aa, alpha in enumerate(alphas):
        plt.plot(x, np.exp(alpha) * x ** 1, memcols[aa], 
                 label = r'$\alpha$ = %4.1f' % alpha)
    plt.xlabel('$x$')
    plt.ylabel('$y$')
    plt.axis([0, 50, 0, 50])
    plt.legend(loc = 4, framealpha = 0.5)
    plt.gca().set_aspect('equal')
    plt.grid()

    plt.subplot(122)
    for bb, beta in enumerate(betas):
        plt.plot(x, np.exp(0) * x ** beta, memcols[bb], 
                 label = r'$\beta$ = %4.2f' % beta)
    plt.xlabel('$x$')
    plt.ylabel('$y$')
    plt.axis([0, 50, 0, 50])
    plt.legend(loc = 1, framealpha = 0.5)
    plt.gca().set_aspect('equal')
    plt.grid()

    plt.savefig('app_mem.pdf')
    plt.close()
