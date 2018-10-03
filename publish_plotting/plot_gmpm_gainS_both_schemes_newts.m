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
sens_ts_figfile_cs = 'GainS_TS_vs_gamma_cs.svg';
sens_ts_figfile_cz = 'GainS_TS_vs_gamma_cz.svg';
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

if 0
[GM_s1, PM_s1, Sens_gain1, TS_s1, BW_s1] = gmpm_vs_gam_recyc_obs(G, G_recyc,...
  G_obsDist, Q_cs, R_cs, S_cs, LxLd, gams);

[GM_s2, PM_s2, Sens_gain2, TS_s2, BW_s2] = gmpm_vs_gam_recyc_obs(G, G_recyc, G_obsDist, Q_cz, R_cz, S_cz, LxLd, gams);

save('gmpmS_data.mat', 'GM_s1', 'PM_s1', 'Sens_gain1', 'TS_s1', 'BW_s1',...
  'GM_s2', 'PM_s2', 'Sens_gain2', 'TS_s2', 'BW_s2') 
else
  load('gmpmS_data.mat')
end
% -------------------------------------------------------------------------
% -------------- Build the LaTeX table ------------------------------------

% ------------------- Constant-sigma --------------------------------------
% 1) Robust-optimal
[~, idx_csro] = min(Sens_gain1);
% for rob-optimal, mpc and linear take the same gamma
lin_sig_rob_data = [gams(idx_csro), GM_s1(idx_csro), PM_s1(idx_csro), Sens_gain1(idx_csro)];
mpc_sig_rob_data = [gams(idx_csro), GM_s1(idx_csro), PM_s1(idx_csro), Sens_gain1(idx_csro)];

% 2). minimum-gamma. For const-sigma, and rmax=14, we get gam=10 for
%     linear, and gam=0.5 for Nmpc=12.
idx_csmg_lin = find(gams == min_gam_slf_cs, 1, 'first');
lin_sig_ming_data = [gams(idx_csmg_lin), GM_s1(idx_csmg_lin),...
                     PM_s1(idx_csmg_lin), Sens_gain1(idx_csmg_lin)];
idx_csmg_mpc = find(gams == min_gam_mpc_cs, 1, 'first');
mpc_sig_ming_data = [gams(idx_csmg_mpc), GM_s1(idx_csmg_mpc),...
                     PM_s1(idx_csmg_mpc), Sens_gain1(idx_csmg_mpc)];

                   
K_cs_rg = dlqr(G_recyc.a, G_recyc.b, Q_cs, R_cs+gams(idx_csro), S_cs);
K_cs_mg_lin = dlqr(G_recyc.a, G_recyc.b, Q_cs, R_cs+gams(idx_csmg_lin), S_cs);
K_cs_mg_mpc = dlqr(G_recyc.a, G_recyc.b, Q_cs, R_cs+gams(idx_csmg_mpc), S_cs);

[Sens_cs_rg]     = ss_loops_delta_dist(G, G_recyc, G_obsDist, K_cs_rg, LxLd);
[Sens_cs_mg_lin] = ss_loops_delta_dist(G, G_recyc, G_obsDist, K_cs_mg_lin, LxLd);
[Sens_cs_mg_mpc] = ss_loops_delta_dist(G, G_recyc, G_obsDist, K_cs_mg_mpc, LxLd);


% --------------------- Choose-zeta ---------------------------------------
% 1). Robust optimal
[~, idx_czro] = min(Sens_gain2);
lin_zet__rob_data = [gams(idx_czro), GM_s2(idx_czro), PM_s2(idx_czro), Sens_gain2(idx_czro)];
mpc_zet_rob_data = [gams(idx_czro), GM_s2(idx_czro), PM_s2(idx_czro), Sens_gain2(idx_czro)];

% 2). We have, I think, gam=0.2 for MPC and gam=2.9 for linear,
%     in the choose zeta scheme.
idx_czmg_lin = find(gams == min_gam_slf_cz, 1, 'first');
lin_zet_ming_data = [gams(idx_czmg_lin), GM_s2(idx_czmg_lin),...
                     PM_s2(idx_czmg_lin), Sens_gain2(idx_czmg_lin)];
idx_czmg_mpc = find(gams == min_gam_mpc_cz, 1, 'first');
mpc_zet_ming_data = [gams(idx_czmg_mpc), GM_s2(idx_czmg_mpc),...
                     PM_s2(idx_czmg_mpc), Sens_gain2(idx_czmg_mpc)];

                   
K_cz_rg = dlqr(G_recyc.a, G_recyc.b, Q_cz, R_cz+gams(idx_czro), S_cz);
K_cz_mg_lin = dlqr(G_recyc.a, G_recyc.b, Q_cz, R_cs+gams(idx_czmg_lin), S_cz);
K_cz_mg_mpc = dlqr(G_recyc.a, G_recyc.b, Q_cz, R_cs+gams(idx_czmg_mpc), S_cz);

[Sens_cz_rg]     = ss_loops_delta_dist(G, G_recyc, G_obsDist, K_cz_rg, LxLd);
[Sens_cz_mg_lin] = ss_loops_delta_dist(G, G_recyc, G_obsDist, K_cz_mg_lin, LxLd);
[Sens_cz_mg_mpc] = ss_loops_delta_dist(G, G_recyc, G_obsDist, K_cz_mg_mpc, LxLd);

%%

clc
width = 3.45;
height = 2.5;
f1 = mkfig(1, width, height); clf
% f1 = figure; clf
freqs = logspace(log10(1), log10(12500), 200);
hsens_cs_rg = frf_bode_mag(Sens_cs_rg, freqs, f1, 'Hz', 'Color', 'k', 'LineStyle', '-');
hsens_cs_rg.DisplayName = 'CS-RG (SLF/MPC)';

hsens_cs_mg_lin = frf_bode_mag(Sens_cs_mg_lin, freqs, f1, 'Hz', 'Color', 'k', 'LineStyle', '--');
hsens_cs_mg_lin.DisplayName = 'CS-MG (SLF)';
hsens_cs_mg_mpc = frf_bode_mag(Sens_cs_mg_mpc, freqs, f1, 'Hz', 'Color', 'k', 'LineStyle', ':');
hsens_cs_mg_mpc.DisplayName = 'CS-MG (MPC)';

hsens_cz_rg = frf_bode_mag(Sens_cz_rg, freqs, f1, 'Hz', 'Color', 'r', 'LineStyle', '-');
hsens_cz_rg.DisplayName = 'CZ-RG (SLF/MPC)';

hsens_cz_mg_lin = frf_bode_mag(Sens_cz_mg_lin, freqs, f1, 'Hz', 'Color', 'r', 'LineStyle', '--');
hsens_cz_mg_lin.DisplayName = 'CZ-MG (SLF)';
hsens_cz_mg_mpc = frf_bode_mag(Sens_cz_mg_mpc, freqs, f1, 'Hz', 'Color', 'r', 'LineStyle', ':');
hsens_cz_mg_mpc.DisplayName = 'CZ-MG (MPC)';

% response of ppure derivitive.
% frf_bode_mag(zpk([1], [], 1, G.Ts), freqs, f1, 'Hz')

leg_sens = legend([hsens_cs_rg, hsens_cs_mg_lin, hsens_cs_mg_mpc,...
  hsens_cz_rg, hsens_cz_mg_lin, hsens_cz_mg_mpc]);

ax = gca();
ax.XTick = [1, 10, 100, 1000, 10000];
ylim(ax, [-50, 15])
tighten_axis(f1, ax)

set(leg_sens, 'FontSize', 7, 'NumColumns', 2, 'Box', 'off');
set(leg_sens, 'Position', [0.5805    0.3410    2.7610    0.4021]);

if saveon
  saveas(f1, fullfile(PATHS.jfig, 'sens_bode.svg'));
end

%
data = {lin_sig_ming_data, 'SLF-CS-MG', 'SLF, const-$\sigma$, minimum-$\gamma$';
  mpc_sig_ming_data, 'MPC-CS-MG', 'MPC, const-$\sigma$ minimum-$\gamma$'; 
  lin_sig_rob_data, 'SLF-CS-RG', 'SLF, const-$\sigma$, robust-$\gamma$ ';
  mpc_sig_rob_data, 'MPC-CS-RG', 'MPC, const-$\sigma$, robust-$\gamma$'; 
  %
  lin_zet_ming_data, 'SLF-CZ-MG', 'SLF, choose-$\zeta$ minimum-$\gamma$'; 
  mpc_zet_ming_data, 'MPC-CZ-MG', 'MPC, choose-$\zeta$ minimum-$\gamma$';
  %
  lin_zet__rob_data, 'SLF-CZ-RG', 'SLF, choose-$\zeta$  robust-$\gamma$'; 
  mpc_zet_rob_data, 'MPC-CZ-RG', 'MPC, choose-$\zeta$ robust-$\gamma$'}


body = sprintf('%s\n\\hline\n', 'scheme    & abbreviation       & $\gamma$ &GM & PM & $|\mathcal{S}|$\\');


for col_cell=data'
  dat = col_cell{1};
  abrev = col_cell{2};
  name = col_cell{3};
  
  row_str = sprintf('%s & %s & ', name, abrev);
  for col = 1:size(dat,2)
    
    if col < size(dat,2)
      fmt = '%s %.1f &';
    else
      fmt = '%s %.1f';
    end
     row_str = sprintf(fmt, row_str, dat(col));
  end
  body = sprintf('%s%s\\\\\n', body, row_str);
end

body

if saveon
  fid = fopen(fullfile(PATHS.MPCJ_root, 'latex', 'rob_data.tex'), 'w+');
  fprintf(fid, '%s', body);
  fclose(fid);
end
%%
clc
% -------------- Plot Colors & Linestyles ----------------------
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
% plot(gams(idx_csro), PM_s1(idx_csro), 'xk')
% 
% yyaxis(ax1, 'left')
% plot(gams(idx_csmg_lin), GM_s1(idx_csmg_lin), 'dk')
% yyaxis(ax1, 'right')
% plot(gams(idx_csmg_lin), PM_s1(idx_csmg_lin), 'dk')
% 
% yyaxis(ax1, 'left')
% plot(gams(idx_csmg_mpc), GM_s1(idx_csmg_mpc), 'sk')
% yyaxis(ax1, 'right')
% plot(gams(idx_csmg_mpc), PM_s1(idx_csmg_mpc), 'sk')

% % choose zeta
% yyaxis(ax1, 'left')
% % plot(gams(idx_czro), GM_s2(idx_czro), 'xk')
% yyaxis(ax1, 'right')
% % plot(gams(idx_czro), PM_s2(idx_czro), 'xk')
% 
% yyaxis(ax1, 'left')
% plot(gams(idx_czmg_lin), GM_s2(idx_czmg_lin), 'dk')
% yyaxis(ax1, 'right')
% plot(gams(idx_czmg_lin), PM_s2(idx_czmg_lin), 'dk')
% 
% yyaxis(ax1, 'left')
%  plot(gams(idx_czmg_mpc), GM_s2(idx_czmg_mpc), 'sk')
% yyaxis(ax1, 'right')
% plot(gams(idx_czmg_mpc), PM_s2(idx_czmg_mpc), 'sk')


leg1 = legend([hgm1, hgm2, hpm1, hpm2]);


if saveon
   saveas(f4, fullfile(PATHS.jfig, pmgm_figfile))
end

%%
clc
% =============================================================
% Load the sum of steps experimental data.
exp_path = fullfile(PATHS.step_exp, 'many_steps_sweep_gamma_21-Sep-2018_01');
ts_total_file = 'ts_totalconst-sig_09-21-2018.mat';

fpath1 = fullfile(exp_path, ts_total_file);
load(fpath1)

gams_ts = ts_sum_exp_results(:,1,1);
ts_sum_exp_results = ts_sum_exp_results(:,2:end,:);
ts_exp_means = mean(ts_sum_exp_results, 3); % average across pages (3rd dim)
ts_exp_stdv = sqrt(var(ts_sum_exp_results,0,3)); % average across pages (3rd dim)
ts_data = struct('means', ts_exp_means, 'stdv', ts_exp_stdv, 'gams', gams_ts);

% ===============================
% figure(f5);
f5 = mkfig(5+figbase, 7, 2.0); clf
ax2 = axes('Units', 'inches', 'Position', [0.3793 0.3000 2.65 1.50],...
  'YAxisLocation', 'left');
[h_sens1] = plot_sens(gams, Sens_gain1, ax2, f5,...
  'Sens_color', GM_colr, 'Sens_LS', GM_ls1);

h_sens1.DisplayName = '$|\mathcal{S}_I(1)|$';

yyaxis left
ax2.YAxis(2).Color = PM_colr;
hts_mpc_exp = errorbar(ts_data.gams, ts_data.means(:,2)*1000,...
  ts_data.stdv(:,2)*1000, 'r.' );
hts_mpc_exp.DisplayName = '$\bar{T}(\gamma_{\ell})$ (MPC exp.)';

hts_slf_exp = errorbar(ts_data.gams, ts_data.means(:,1)*1000,...
  ts_data.stdv(:,1)*1000, 'k.' );
hts_slf_exp.DisplayName = '$\bar{T}(\gamma_{\ell})$ (SLF exp.)';

title('constant-$\sigma$')

hold on
ylabel('total settle-time [ms]')
yyaxis(ax2, 'left')
ylim([200, 300])

xlim(xlm)
yyaxis(ax2, 'right')

set(ax2, 'XTick', [0.0001, 0.001, 0.01, 0.1, 1, 10, 100]);
% set(ax2, 'XTickLabel', []);
xlab1 = xlabel('$\gamma$');
set(xlab1, 'Units', 'inches', 'Position', [1.349, -.15, 0])
set(ax2, 'YTick', [36:1:39]);
% ----------------------------------------------------------------------- %
exp_path = fullfile(PATHS.step_exp, 'many_steps_sweep_gamma_21-Sep-2018_01');
ts_total_file = 'ts_totalchoose-zet_09-21-2018.mat';

fpath2 = fullfile(exp_path, ts_total_file);
load(fpath2)

gams_ts = ts_sum_exp_results(:,1,1);

ts_sum_exp_results = ts_sum_exp_results(:,2:end,:);
ts_exp_means = mean(ts_sum_exp_results, 3); % average across pages (3rd dim)
ts_exp_stdv = sqrt(var(ts_sum_exp_results,0,3)); % average across pages (3rd dim)
ts_data = struct('means', ts_exp_means, 'stdv', ts_exp_stdv, 'gams', gams_ts);

% ===============================
% figure(f6);

ax3 = axes('Units', 'inches', 'Position', [0.3793+3.54 0.3000 2.65 1.50],...
  'YAxisLocation', 'left');
[h_sens2] = plot_sens(gams, Sens_gain2, ax3, f5,...
  'Sens_color', GM_colr, 'Sens_LS', GM_ls1);
h_sens2.DisplayName = '$|\mathcal{S}_I(1)|$';

yyaxis left
ax3.YAxis(2).Color = PM_colr;
hts_mpc_exp2 = errorbar(ts_data.gams, ts_data.means(:,2)*1000,...
  ts_data.stdv(:,2)*1000, 'r.' );
hts_mpc_exp2.DisplayName = '$\bar{T}(\gamma_{\ell})$ (MPC exp.)';

hts_slf_exp2 = errorbar(ts_data.gams, ts_data.means(:,1)*1000,...
  ts_data.stdv(:,1)*1000, 'k.' );
hts_slf_exp2.DisplayName = '$\bar{T}(\gamma_{\ell})$ (SLF exp.)';

hold on
% h_ts_slf_sim2 = plot(gams_ts, ts_sum_sim_results(:,2)*1000, 'k.');
% h_ts_slf_sim2.DisplayName = '$\Sigma t_{s}^i$ (sim.)';
% 
% h_ts_mpc_sim2 = plot(gams_ts, ts_sum_sim_results(:,3)*1000, 'ro', 'MarkerSize', 5);
% h_ts_mpc_sim2.DisplayName = '$\Sigma t_{s}^i$ (sim.)';
ylabel('total settle-time [ms]')

title('choose-$\zeta$')
yyaxis(ax3, 'right')
xlim(xlm)
ylim([44, 48.2])

yyaxis(ax3, 'left')
ylim([200, 700])

leg22 = legend([hts_mpc_exp, hts_slf_exp, h_sens1]);
leg22.Position =  [0.0649 0.6289 0.1875 0.2285];

set(ax3, 'XTick', [0.0001, 0.001, 0.01, 0.1, 1, 10, 100])

xlab2 = xlabel('$\gamma$');
set(xlab2, 'Units', 'inches', 'Position', [1.349, -.15, 0])
if saveon
  saveas(f5, fullfile(PATHS.jfig, sens_ts_figfile_both ))
end

%%


dat1 = load(fpath1);
dat2 = load(fpath2);
f10 = mkfig(10+figbase, 7, 2.0); clf
ax2 = axes('Units', 'inches', 'Position', [0.3793 0.3000 2.65 1.50],...
  'YAxisLocation', 'left');


h_ts_slf_sim = semilogx(dat1.gam_s, dat1.ts_sum_sim_results(:,2)*1000, 'k.');
h_ts_slf_sim.DisplayName = '$\Sigma t_i (\gamma_{\ell})$  (SLF sim.)';
hold on
h_ts_mpc_sim = semilogx(dat1.gam_s, dat1.ts_sum_sim_results(:,3)*1000, 'ro');
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
set(ax2, 'XTick', [0.0001, 0.001, 0.01, 0.1, 1, 10, 100])
grid on

% subplot(1,2,2)
ax3 = axes('Units', 'inches', 'Position', [0.3793+3.54 0.3000 2.65 1.50],...
  'YAxisLocation', 'left');

h_ts_slf_sim = semilogx(dat2.gam_s, dat2.ts_sum_sim_results(:,2)*1000, 'k.');
h_ts_slf_sim.DisplayName = '$\Sigma t_i(\gamma_{\ell})$  (SLF sim.)';
hold on
h_ts_mpc_sim = semilogx(dat2.gam_s, dat2.ts_sum_sim_results(:,3)*1000, 'ro');
h_ts_mpc_sim.DisplayName = '$\Sigma t_i(\gamma_{\ell})$ (MPC sim.)';
ylabel('total settle-time [ms]')
ylim([100, 600])
yyaxis right
ax3.YAxis(2).Color = 'k';
h_bw = semilogx(gams, BW_s2/2/pi, '-k');
h_bw.DisplayName = 'Closed-loop Bandwidth';

title('choose-$\zeta$')
ylabel('Bandwidth [Hz]')
xlab2 = xlabel('$\gamma$');
set(xlab2, 'Units', 'inches', 'Position', [1.349, -.15, 0])
set(ax3, 'XTick', [0.0001, 0.001, 0.01, 0.1, 1, 10, 100])

leg3 = legend([h_ts_mpc_sim, h_ts_slf_sim, h_bw]);
leg3.Position = [0.5675 0.4467 0.2405 0.2285];

xlim(xlm)
grid on
if saveon
  saveas(f10, fullfile(PATHS.jfig, bw_ts_figfile_both ))
end
%%
% clc
[cor_cs_sim_mpc, p_cs_sim_mpc] = corrcoef(ts_sum_sim_results(:,3), ts_data.gams)
[cor_cs_exp_mpc, p_cs_exp_mpc] = corrcoef(ts_data.means(1:end-2,2), ts_data.gams(1:end-2))

sg2 = 10.^(Sens_gain2/20);

sens_interp = interp1(gams, sg2, ts_data.gams);

[cor_cs_sim_mpc, p_cs_sim_mpc] = corrcoef(ts_data.means(:,2), sens_interp)

figure(104)
plot(sens_interp, ts_data.means(:,2), 'x')


figure(105), clf
hold on, grid on
plot(ts_data.gams, ts_data.means(:,2))
yyaxis right
plot(ts_data.gams, sens_interp)


%%
% % % %%
% % % clc
% % % % =============================================================
% % % % Load the sum of steps experimental data.
% % % % exp_path = fullfile(PATHS.step_exp, 'many_steps_sweep_gamma_17-Sep-2018_01');
% % % % ts_total_file = 'ts_totalconst-sig_09-17-2018.mat';
% % % exp_path = fullfile(PATHS.step_exp, 'many_steps_sweep_gamma_21-Sep-2018_01');
% % % ts_total_file = 'ts_totalconst-sig_09-21-2018.mat';
% % % 
% % % 
% % % 
% % % fpath = fullfile(exp_path, ts_total_file);
% % % load(fpath)
% % % 
% % % gams_ts = ts_sum_exp_results(:,1,1);
% % % ts_sum_exp_results = ts_sum_exp_results(:,2:end,:);
% % % ts_exp_means = mean(ts_sum_exp_results, 3); % average across pages (3rd dim)
% % % ts_exp_stdv = sqrt(var(ts_sum_exp_results,0,3)); % average across pages (3rd dim)
% % % ts_data = struct('means', ts_exp_means, 'stdv', ts_exp_stdv, 'gams', gams_ts);
% % % 
% % % 
% % % % ===============================
% % % % figure(f5);
% % % f5 = mkfig(5+figbase, width, 3.66); clf
% % % % [0.1100 0.1300 0.750 0.690]
% % % ax2 = axes('Position', [0.1100 0.4850 0.7800 0.3741], 'YAxisLocation', 'left');
% % % [h_sens1] = plot_sens(gams, Sens_gain1, ax2, f5,...
% % %   'Sens_color', GM_colr, 'Sens_LS', GM_ls1);
% % % % set(ax2, 'XTick', [0.01, 0.1, 1, 10, 100])
% % % % % [h_sens1, h_ts1] = plot_sens_ts(gams, Sens_gain1, TS_s1, ax2, f5,...
% % % % %   'Ts_color', PM_colr, 'Ts_LS', PM_ls1, 'Sens_color', GM_colr, 'Sens_LS', GM_ls1);
% % % 
% % % h_sens1.DisplayName = '$|\mathcal{S}_I(1)|$';
% % % 
% % % yyaxis right
% % % ax2.YAxis(2).Color = PM_colr;
% % % hts_mpc_exp = errorbar(ts_data.gams, ts_data.means(:,2)*1000,...
% % %   ts_data.stdv(:,2)*1000, 'r.' );
% % % hts_mpc_exp.DisplayName = '$\Sigma t_{s}^i$ (MPC exp.)';
% % % 
% % % hts_slf_exp = errorbar(ts_data.gams, ts_data.means(:,1)*1000,...
% % %   ts_data.stdv(:,1)*1000, 'k.' );
% % % hts_slf_exp.DisplayName = '$\Sigma t_{s}^i$ (SLF exp.)';
% % % 
% % % 
% % % hold on
% % % h_ts_slf_sim = plot(gams_ts, ts_sum_sim_results(:,2)*1000, 'k.');
% % % h_ts_slf_sim.DisplayName = '$\Sigma t_{s}^i$  (SLF sim.)';
% % % 
% % % h_ts_mpc_sim = plot(gams_ts, ts_sum_sim_results(:,3)*1000, 'ro');
% % % h_ts_mpc_sim.DisplayName = '$\Sigma t_{s}^i$ (MPC sim.)';
% % % ylabel('total settle-time [ms]')
% % % yyaxis(ax2, 'right')
% % % ylim([100, 275])
% % % % xlim([0.0001, 1000])
% % % xlim(xlm)
% % % yyaxis(ax2, 'left')
% % % % ylim([36, 48])
% % % 
% % % % plot(gams(idx_csro), Sens_gain1(idx_csro), 'xk');
% % % plot(gams(idx_csmg_lin), Sens_gain1(idx_csmg_lin), 'dk');
% % % plot(gams(idx_csmg_mpc), Sens_gain1(idx_csmg_mpc), 'sk');
% % % 
% % % % leg2 = legend([h_sens1, hts_mpc_exp, h_ts_sim]);
% % % % set(leg2, 'Position', [0.1514 0.8066 0.7184 0.2091], 'Box', 'off')
% % % set(ax2, 'XTick', [0.0001, 0.001, 0.01, 0.1, 1, 10, 100]);
% % % set(ax2, 'XTickLabel', []);
% % % xlabel('');
% % % set(ax2, 'YTick', [36:1:39]);
% % % 
% % % 
% % % % exp_path = fullfile(PATHS.step_exp, 'many_steps_sweep_gamma_17-Sep-2018_01');
% % % % ts_total_file = 'ts_totalchoose-zet_09-17-2018.mat';
% % % exp_path = fullfile(PATHS.step_exp, 'many_steps_sweep_gamma_21-Sep-2018_01');
% % % ts_total_file = 'ts_totalchoose-zet_09-21-2018.mat';
% % % 
% % % fpath = fullfile(exp_path, ts_total_file);
% % % load(fpath)
% % % 
% % % gams_ts = ts_sum_exp_results(:,1,1);
% % % 
% % % ts_sum_exp_results = ts_sum_exp_results(:,2:end,:);
% % % ts_exp_means = mean(ts_sum_exp_results, 3); % average across pages (3rd dim)
% % % ts_exp_stdv = sqrt(var(ts_sum_exp_results,0,3)); % average across pages (3rd dim)
% % % ts_data = struct('means', ts_exp_means, 'stdv', ts_exp_stdv, 'gams', gams_ts);
% % % 
% % % % ===============================
% % % % figure(f6);
% % % % f6 = mkfig(6+figbase, width, height); clf
% % % ax3 = axes('Position', [0.1100 0.0856 0.7800 0.3619], 'YAxisLocation', 'left');
% % % [h_sens2] = plot_sens(gams, Sens_gain2, ax3, f5,...
% % %   'Sens_color', GM_colr, 'Sens_LS', GM_ls1);
% % % h_sens2.DisplayName = '$|\mathcal{S}_I(1)|$';
% % % % plot(gams(idx_czro), Sens_gain2(idx_czro), 'xk')
% % % plot(gams(idx_czmg_lin), Sens_gain2(idx_czmg_lin), 'dk')
% % % plot(gams(idx_czmg_mpc), Sens_gain2(idx_czmg_mpc), 'sk')
% % % 
% % % 
% % % yyaxis right
% % % ax3.YAxis(2).Color = PM_colr;
% % % hts_mpc_exp2 = errorbar(ts_data.gams, ts_data.means(:,2)*1000,...
% % %   ts_data.stdv(:,2)*1000, 'r.' );
% % % hts_mpc_exp2.DisplayName = '$\Sigma t_{s}^i$ (MPC exp.)';
% % % 
% % % hts_slf_exp2 = errorbar(ts_data.gams, ts_data.means(:,1)*1000,...
% % %   ts_data.stdv(:,1)*1000, 'k.' );
% % % hts_slf_exp2.DisplayName = '$\Sigma t_{s}^i$ (SLF exp.)';
% % % 
% % % % yyaxis right
% % % % 
% % % % ax2.YAxis(2).Color = PM_colr;
% % % % h_ts_exp2 = errorbar(ts_data.gams, ts_data.means(:,2)*1000,...
% % % %   ts_data.stdv(:,2)*1000, 'r.' );
% % % % h_ts_exp2.DisplayName = 'total settle-time (experiment)';
% % % hold on
% % % h_ts_slf_sim2 = plot(gams_ts, ts_sum_sim_results(:,2)*1000, 'k.');
% % % h_ts_slf_sim2.DisplayName = '$\Sigma t_{s}^i$ (sim.)';
% % % 
% % % h_ts_mpc_sim2 = plot(gams_ts, ts_sum_sim_results(:,3)*1000, 'ro', 'MarkerSize', 5);
% % % h_ts_mpc_sim2.DisplayName = '$\Sigma t_{s}^i$ (sim.)';
% % % ylabel('total settle-time [ms]')
% % % 
% % % 
% % % yyaxis(ax3, 'left')
% % % % ylim([36, 48])
% % % xlim(xlm)
% % % ylim([42, 48.2])
% % % 
% % % 
% % % yyaxis(ax3, 'right')
% % % ylim([00, 700])
% % % 
% % % leg22 = legend([h_sens1, hts_mpc_exp, hts_slf_exp, h_ts_mpc_sim, h_ts_slf_sim]);
% % % set(leg22, 'NumColumns', 2, 'Box', 'off', 'Position', [0.1189 0.8784 0.7491 0.1152]);
% % % 
% % % 
% % % 
% % % set(ax3, 'XTick', [0.0001, 0.001, 0.01, 0.1, 1, 10, 100])
% % % % set(ax2, 'YTick', [36:1:39])
% % % 
% % % 
% % % % legend([h_sens1, h_sens2, h_ts1, h_ts2]);
%%



if 1
   saveas(f4, fullfile(PATHS.jfig, pmgm_figfile))
   saveas(f5, fullfile(PATHS.jfig, sens_ts_figfile_cs))
%    saveas(f6, fullfile(PATHS.jfig, sens_ts_figfile_cz))
end

%



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




