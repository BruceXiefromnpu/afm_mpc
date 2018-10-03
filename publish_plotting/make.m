% This script acts like a "makefile" to build all the plots.

% Figure 3, frequency response and fit for stage in x-direction
plot_FRF_xdir

% Table I (state-space data with small pz-plot)
table_pzdata

plot_fit_drift
plot_hyst_sat_fit
%%
plot_lqr_locus_constsigma
plot_lqr_locus_choosezeta
