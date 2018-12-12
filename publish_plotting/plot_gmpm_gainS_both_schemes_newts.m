% This script plots the LQR based root locus for the "lowgain"
% version of things.

addpath(fullfile(getMatPath(), 'afm_mpc_journal', 'functions'))
addpath(fullfile(getMatPath(), 'afm_mpc_journal', 'functions', 'canon'))
addpath('~/gradschool/sysID/matlab/functions/')

clear, clc
saveon = true;
plants = CanonPlants.plants_ns14(9,2);
G = plants.SYS;
G_recyc = plants.sys_recyc;
Ts = G_recyc.Ts;

pmgm_figfile = 'PMGM_vs_gamma_both.svg';
sens_ts_figfile_both = 'GainS_TS_vs_gamma_both.svg';
bw_ts_figfile_both = 'BW_TS_vs_gamma_both.svg';

Qw = plants.SYS.b*plants.SYS.b'*50;
Lx = G.a*dlqr(plants.SYS.a', plants.SYS.c', Qw, 1)';
p_int_d = 0.7;
[LxLd, G_obsDist, Ident_obs, C_ydist] = DistEst.output_dist_est(G, Lx, p_int_d);

cmplx_rad = 0.9;
% Constant sigma LQR weights
[Q_cs, R_cs, S_cs, P_x] = build_control_constsigma(G_recyc, cmplx_rad);
% Chooze zeta LQR weights
can_cntrl = CanonCntrlParams_ns14();
[Q_cz, R_cz, S_cz] = build_control(G_recyc, can_cntrl);

gams = logspace(log10(0.1), log10(1000), 100);
% gams = unique([gams, 1.0, 3.0, 3.1, 3.8, 12,10, 2, 0.001, 0.5, 0.2, 2.9]);
min_gam_slf_cs = 7.5;
min_gam_mpc_cs = 0.001;
min_gam_slf_cz = 3.5;
min_gam_mpc_cz = 0.00001;
gams = unique([gams, 1.0, 3.0, 3.1, 3.8, min_gam_slf_cs, min_gam_mpc_cs,...
  min_gam_slf_cz, min_gam_mpc_cz]);

if 1
  [GM_s1, PM_s1, Sens_gain1, TS_s1, BW_s1] = gmpm_vs_gam_recyc_obs(G, G_recyc,...
    G_obsDist, Q_cs, R_cs, S_cs, LxLd, gams);
  
  [GM_s2, PM_s2, Sens_gain2, TS_s2, BW_s2] = gmpm_vs_gam_recyc_obs(G, G_recyc,...
    G_obsDist, Q_cz, R_cz, S_cz, LxLd, gams);
  
%   save('gmpmS_data.mat', 'GM_s1', 'PM_s1', 'Sens_gain1', 'TS_s1', 'BW_s1',...
%     'GM_s2', 'PM_s2', 'Sens_gain2', 'TS_s2', 'BW_s2')
else
  load('gmpmS_data.mat')
end


% --------------- The PM-GM figure (fig 12) --------------------

%  Plot Colors & Linestyles ----------------------

GM_colr = 'k';
PM_colr = 'r';
GM_ls1 = '-';
GM_ls2 = '--';

PM_ls1 = '-.';
PM_ls2 = ':';

width = 3.45;
height = 2.5;
figbase = 0;
f4 = mkfig(4+figbase, width, height); clf
figure(f4);
ax1 = axes('Position', [0.1100 0.1300 0.7750 0.8150], 'YAxisLocation', 'left');

[hgm1, hpm1] = plot_gmpm(gams, GM_s1, PM_s1, ax1, f4, 'GM_color', GM_colr,...
  'GM_LS', GM_ls1, 'PM_color', PM_colr, 'PM_LS', PM_ls1);
hgm1.DisplayName = 'GM (constant-$\sigma$)';
hpm1.DisplayName = 'PM (constant-$\sigma$)';

xlm = [0.000005, 1000]; % x-limits for all figs
figure(f4);
[hgm2, hpm2] = plot_gmpm(gams, GM_s2, PM_s2, ax1, f4, 'GM_color', GM_colr,...
  'GM_LS', GM_ls2, 'PM_color', PM_colr, 'PM_LS', PM_ls2);
hgm2.DisplayName = 'GM (choose-$\zeta$)';
hpm2.DisplayName = 'PM (choose-$\zeta$)';

xlim(ax1, xlm);
set(ax1, 'XTick', [.0001, .001, 0.01, 0.1, 1, 10, 100])
leg1 = legend([hgm1, hgm2, hpm1, hpm2]);
set(leg1, 'Position', [0.1374 0.6985 0.4108 0.2509], 'Box', 'off')
%

f4.CurrentAxes = ax1;
yyaxis(ax1, 'left')
% plot(gams(idx_csro), GM_s1(idx_csro), 'xk')
yyaxis(ax1, 'right')


leg1 = legend([hgm1, hgm2, hpm1, hpm2]);


if saveon
   saveas(f4, fullfile(PATHS.jfig, pmgm_figfile))
end

%
% ----------------------------------------------------------------------- %
% Now, plot the settling settling times vs gamma and on the same two graphs,
% plot the DC-gain of the integrated sensitivity function vs gamma.
% ----------------------------------------------------------------------- %

%%
% =============================================================
% Load the CONSTANT-SIMGA sum of steps experimental data.
exp_path_CS = fullfile(PATHS.step_exp, 'many_steps_sweep_gamma_21-Sep-2018_01');
ts_total_file_CS = 'ts_totalconst-sig_09-21-2018.mat';

fpath1 = fullfile(exp_path_CS, ts_total_file_CS);
dat_CS=load(fpath1);

% gams_ts = ts_sum_exp_results(:,1,1);
ts_sum_exp_results_CS = dat_CS.ts_sum_exp_results(:,2:end,:); % first column is gammas, drop it.
ts_exp_means_CS = mean(ts_sum_exp_results_CS, 3); % average across pages (3rd dim)
ts_exp_stdv_CS = sqrt(var(ts_sum_exp_results_CS,0,3)); % average across pages (3rd dim)
ts_data_CS = struct('means', ts_exp_means_CS, 'stdv', ts_exp_stdv_CS, 'gams', dat_CS.gam_s);

% Load the CHOOSE-ZETA sum of steps experimental data.
exp_path_CZ = fullfile(PATHS.step_exp, 'many_steps_sweep_gamma_21-Sep-2018_01');
ts_total_file_CZ = 'ts_totalchoose-zet_09-21-2018.mat';
fpath2 = fullfile(exp_path_CZ, ts_total_file_CZ);
dat_CZ = load(fpath2);

%gams_ts = ts_sum_exp_results(:,1,1);
ts_sum_exp_results_CZ = dat_CZ.ts_sum_exp_results(:,2:end,:);
ts_exp_means_CZ = mean(ts_sum_exp_results_CZ, 3); % average across pages (3rd dim)
ts_exp_stdv_CZ = sqrt(var(ts_sum_exp_results_CZ,0,3)); % average across pages (3rd dim)
ts_data_CZ = struct('means', ts_exp_means_CZ, 'stdv', ts_exp_stdv_CZ, 'gams', dat_CZ.gam_s);

% =============================================================

f5 = mkfig(5+figbase, 7, 2.0); clf
ax2 = axes('Units', 'inches', 'Position', [0.3793 0.3000 2.65 1.50],...
  'YAxisLocation', 'left');

ax3 = axes('Units', 'inches', 'Position', [0.3793+3.54 0.3000 2.65 1.50],...
  'YAxisLocation', 'left');

% --------------------- plot CONSTANT-SIGMA ----------------------------- %
f5.CurrentAxes = ax2;
[h_sens1] = plot_sens(gams, Sens_gain1, ax2, f5,...
  'Sens_color', GM_colr, 'Sens_LS', GM_ls1);
h_sens1.DisplayName = '$|\mathcal{S}_I(1)|$';

yyaxis left
ax2.YAxis(2).Color = PM_colr;
hts_mpc_exp = errorbar(ts_data_CS.gams, ts_data_CS.means(:,2)*1000,...
  ts_data_CS.stdv(:,2)*1000, 'r.' );
hts_mpc_exp.DisplayName = '$\bar{T}(\gamma_{\ell})$ (MPC exp.)';

hts_slf_exp = errorbar(ts_data_CS.gams, ts_data_CS.means(:,1)*1000,...
  ts_data_CS.stdv(:,1)*1000, 'k.' );
hts_slf_exp.DisplayName = '$\bar{T}(\gamma_{\ell})$ (SLF exp.)';

title('constant-$\sigma$')

hold on
ylabel('total settle-time [ms]')
yyaxis(ax2, 'left')
ylim([200, 300])

xlim(xlm)
yyaxis(ax2, 'right')
ax2.YAxis(2).Color = 'k';

set(ax2, 'XTick', [0.00001, 0.0001, 0.001, 0.01, 0.1, 1, 10, 100]);
xlab1 = xlabel('$\gamma$');
set(xlab1, 'Units', 'inches', 'Position', [1.349, -.15, 0])
set(ax2, 'YTick', [36:1:39]);

% --------------------- plot CHOOSE-ZETA -------------------------------- %
f5.CurrentAxes = ax3;
[h_sens2] = plot_sens(gams, Sens_gain2, ax3, f5,...
  'Sens_color', GM_colr, 'Sens_LS', GM_ls1);
h_sens2.DisplayName = '$|\mathcal{S}_I(1)|$';

yyaxis left
ax3.YAxis(2).Color = PM_colr;
hts_mpc_exp2 = errorbar(ts_data_CZ.gams, ts_data_CZ.means(:,2)*1000,...
  ts_data_CZ.stdv(:,2)*1000, 'r.' );
hts_mpc_exp2.DisplayName = '$\bar{T}(\gamma_{\ell})$ (MPC exp.)';

hts_slf_exp2 = errorbar(ts_data_CZ.gams, ts_data_CZ.means(:,1)*1000,...
  ts_data_CZ.stdv(:,1)*1000, 'k.' );
hts_slf_exp2.DisplayName = '$\bar{T}(\gamma_{\ell})$ (SLF exp.)';

hold on

ylabel('total settle-time [ms]')

title('choose-$\zeta$')
yyaxis(ax3, 'right')
ax3.YAxis(2).Color = 'k';
xlim(xlm)
ylim([44, 48.2])

yyaxis(ax3, 'left')
ylim([200, 700])

leg22 = legend([hts_mpc_exp, hts_slf_exp, h_sens1]);
leg22.Position =  [0.0649 0.6289 0.1875 0.2285];

set(ax3, 'XTick', [1e-5, 0.0001, 0.001, 0.01, 0.1, 1, 10, 100])

xlab2 = xlabel('$\gamma$');
set(xlab2, 'Units', 'inches', 'Position', [1.349, -.15, 0])
if saveon
  saveas(f5, fullfile(PATHS.jfig, sens_ts_figfile_both ))
end

%
% ----------------------------------------------------------------------- %
% Now, plot the simulated settling times vs gamma and on the same two graphs,
% plot the 3-dB closed loop bandwidth vs gamma
% ----------------------------------------------------------------------- %

f10 = mkfig(10+figbase, 7, 2.0); clf
ax2 = axes('Units', 'inches', 'Position', [0.3793 0.3000 2.65 1.50],...
  'YAxisLocation', 'left');


h_ts_slf_sim = semilogx(dat_CS.gam_s, dat_CS.ts_sum_sim_results(:,2)*1000, 'k.');
h_ts_slf_sim.DisplayName = '$\Sigma t_i (\gamma_{\ell})$  (SLF sim.)';
hold on
h_ts_mpc_sim = semilogx(dat_CS.gam_s, dat_CS.ts_sum_sim_results(:,3)*1000, 'ro');
h_ts_mpc_sim.DisplayName = '$\Sigma t_i(\gamma_{\ell})$ (MPC sim.)';
ylabel('total settle-time [ms]')
xlab1 = xlabel('$\gamma$');
set(xlab1, 'Units', 'inches', 'Position', [1.349, -.15, 0])

yyaxis right
ax2.YAxis(2).Color = 'k';
semilogx(gams, BW_s1/2/pi, '-k')
title('constant-$\sigma$')
ylabel('Bandwidth [Hz]')
xlim(xlm)
set(ax2, 'XTick', [1e-5, 0.0001, 0.001, 0.01, 0.1, 1, 10, 100])
grid on

% subplot(1,2,2)
ax3 = axes('Units', 'inches', 'Position', [0.3793+3.54 0.3000 2.65 1.50],...
  'YAxisLocation', 'left');

h_ts_slf_sim = semilogx(dat_CZ.gam_s, dat_CZ.ts_sum_sim_results(:,2)*1000, 'k.');
h_ts_slf_sim.DisplayName = '$\Sigma t_i(\gamma_{\ell})$  (SLF sim.)';
hold on
h_ts_mpc_sim = semilogx(dat_CZ.gam_s, dat_CZ.ts_sum_sim_results(:,3)*1000, 'ro');
h_ts_mpc_sim.DisplayName = '$\Sigma t_i(\gamma_{\ell})$ (MPC sim.)';
ylabel('total settle-time [ms]')
ylim([100, 650])
yyaxis right
ax3.YAxis(2).Color = 'k';
h_bw = semilogx(gams, BW_s2/2/pi, '-k');
h_bw.DisplayName = 'Closed-loop Bandwidth';

title('choose-$\zeta$')
ylabel('Bandwidth [Hz]')
xlab2 = xlabel('$\gamma$');
set(xlab2, 'Units', 'inches', 'Position', [1.349, -.15, 0])
set(ax3, 'XTick', [1e-5, 0.0001, 0.001, 0.01, 0.1, 1, 10, 100])

leg3 = legend([h_ts_mpc_sim, h_ts_slf_sim, h_bw]);
leg3.Position = [0.5675 0.4467 0.2405 0.2285];

xlim(xlm)
grid on
if saveon
  saveas(f10, fullfile(PATHS.jfig, bw_ts_figfile_both ))
end
%%


% ------------------------------------------------------------------
% ------------------------------------------------------------------
% ------- Function Defs. So nice matlab FINALLY allows this --------

function [h_sens] = plot_sens(gams, Sens_gain_s, ax, Fig, varargin)
  p = inputParser();
  p.addParameter('Sens_color', 'k')
  p.addParameter('Sens_LS', '-')
  p.parse(varargin{:})
  Sens_opts = {'Color', p.Results.Sens_color, 'LineStyle', p.Results.Sens_LS};

  
  
  Fig.CurrentAxes = ax;
  yyaxis right;
  ax.YAxis(1).Color = p.Results.Sens_color;
  
  h_sens = semilogx(gams, Sens_gain_s, Sens_opts{:});
  %   h_sens = plot(gams, Sens_gain_s, Sens_opts{:});
  hold on
  grid on
  xlabel('$\gamma$')
  ylabel('$|\mathcal{S}_I(1)|$ [dB]')
  
%   yyaxis right
%   h_ts = semilogx(ts_data.gams, ts_data.means(:,2)*1000, Ts_opts{:});
%   ylabel('Settle-time [ms]')
%   
%   ax.YAxis(2).Color = p.Results.Ts_color;
  
end




function [h_sens, h_ts] = plot_sens_ts(gams, Sens_gain_s, TS_s, ax, Fig, varargin)
  p = inputParser();
  p.addParameter('Sens_color', 'k')
  p.addParameter('Ts_color', 'r')
  p.addParameter('Sens_LS', '-')
  p.addParameter('Ts_LS', '--')
  p.parse(varargin{:})
  Sens_opts = {'Color', p.Results.Sens_color, 'LineStyle', p.Results.Sens_LS};
  Ts_opts = {'Color', p.Results.Ts_color, 'LineStyle', p.Results.Ts_LS};
  ax.YAxis(1).Color = p.Results.Sens_color;
  
  Fig.CurrentAxes = ax;
  yyaxis left;
  
  h_sens = semilogx(gams, Sens_gain_s, Sens_opts{:});
  hold on
  grid on
  xlabel('$\gamma$')
  ylabel('$|\mathcal{S}(z)/(z-1)|_{z=1}|$ [dB]')
  
  yyaxis right
  h_ts = semilogx(gams, TS_s*1000, Ts_opts{:});
  ylabel('Settle-time [ms]')
  
  ax.YAxis(2).Color = p.Results.Ts_color;
  
end


function [h_gm, h_pm] = plot_gmpm(gams, GM_s, PM_s, ax1, Fig, varargin)
  %PLOT_GMPM Summary of this function goes here
  %   Detailed explanation goes here

  p = inputParser();
  p.addParameter('GM_color', 'k')
  p.addParameter('PM_color', 'r')
  p.addParameter('GM_LS', '-')
  p.addParameter('PM_LS', '--')
  p.parse(varargin{:})
  GM_opts = {'Color', p.Results.GM_color, 'LineStyle', p.Results.GM_LS};
  PM_opts = {'Color', p.Results.PM_color, 'LineStyle', p.Results.PM_LS};
  
  
  Fig.CurrentAxes = ax1;
  yyaxis left;
  h_gm = semilogx(ax1, gams, GM_s, GM_opts{:});
  h_gm.DisplayName = 'Gain Margin';
  hold on
  grid on
  xlab = xlabel('$\gamma$');
  ylabel('Gain Margin [dB]')
  ax1.YAxis(1).Color = p.Results.GM_color;
  
  Fig.CurrentAxes = ax1;
  yyaxis right
  h_pm = plot(ax1, gams, PM_s, PM_opts{:});
  h_pm.DisplayName = 'Phase Margin';
  ylabel('Phase Margin [deg]')
  
  set(xlab, 'Units', 'normalized', 'Position', [0.5, -0.051, 0],...
    'FontSize', 14);
  

  xlim([gams(1), gams(end)])
  ax1.YAxis(2).Color = p.Results.PM_color;
  
end




