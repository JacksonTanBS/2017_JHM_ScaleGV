The text files within provide evaluating variables as a function of spatial scale, temporal scale, and threshold. The files are named according to:

[dataset]_[variable]_[threshold].txt

The variables are:

 1) hits (hit)
 2) misses (miss)
 3) false alarms (false)
 4) correct negatives (corneg)
 5) probability of detection (pod)
 6) false alarm ratio (far)
 7) bias in detection (bid)
 8) Heidke skill score (hss)
 9) correlation coefficient (corr)
10) normalized mean error (nme)
11) normalized mean absolute error (nmae)
12) normalized root-mean-square error (nrmse)
13) alpha from the multiplicative error model (mem-alpha)
14) beta from the multiplicative error model (mem-beta)
15) sigma from the multiplicative error model (mem-sigma)

The threshold range from 0 mm / h to 1 mm / h at 0.05 mm / h increments. Note that the 0 mm / h threshold is actually calculated with 1e-6 mm / h, so zero rain rates are still considered not raining.

Each file correspond to the threshold it is calculated with. Within each file is a table. Each column represents increasing temporal resolution: for IMERG, they correspond to 0.5h, 1 h, 3 h, 6 h, 12 h and 24 h; for TMPA, they correspond to 3 h, 6 h, 12 h and 24 h. Each row represents increasing spatial resolution: for IMERG, they correspond to 0.1°, 0.2°, ..., 2.5°; for TMPA, they correspond to 0.25°, 0.50°, ..., 2.50°.

Jackson Tan
jackson.tan@nasa.gov
Oct 31, 2016
