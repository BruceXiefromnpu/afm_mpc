% This script plots the LQR based root locus for the "lowgain"
% version of things.

addpath(fullfile(getMatPath(), 'afm_mpc_journal', 'functions'))
addpath(fullfile(getMatPath(), 'afm_mpc_journal', 'functions', 'canon'))

clear, clc
saveon = true;
plants = CanonPlants.plants_ns14(9, 2);
G = plants.SYS;
G_recyc = plants.sys_recyc;
Ts = G_recyc.Ts;


pmgm_figfile = 'PMGM_vs_gamma_both.svg';
sens_ts_figfile = 'GainS_TS_vs_gamma_both.svg';
lqr_locus_figfile = 'lqr_locus_both.svg';

const_sig = false;


Qw = plants.SYS.b*plants.SYS.b'*50;
Lx = G.a*dlqr(plants.SYS.a', plants.SYS.c', Qw, 1)';
p_int_d = 0.7;
[LxLd, G_obsDist, Ident_obs, C_ydist] = DistEst.output_dist_est(G, Lx, p_int_d);

cmplx_rad = 0.9;
[Q_cs, R_cs, S_cs, P_x] = build_control_constsigma(G_recyc, cmplx_rad);
   

gams = logspace(log10(0.1), log10(1000), 100);
gams = unique([gams, 1.0, 3.0, 3.1, 3.8, 12,10, 2, 0.001, 0.5, 0.2, 2.9]);
[GM_s1, PM_s1, Sens_gain1, TS_s1] = gmpm_vs_gam_recyc_obs(G, G_recyc,...
  G_obsDist, Q_cs, R_cs, S_cs, LxLd, gams);

can_cntrl = CanonCntrlParams_ns14();
[Q_cz, R_cz, S_cz] = build_control(G_recyc, can_cntrl);
[GM_s2, PM_s2, Sens_gain2, TS_s2] = gmpm_vs_gam_recyc_obs(G, G_recyc, G_obsDist, Q_cz, R_cz, S_cz, LxLd, gams);





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
idx_csmg_lin = find(gams == 10, 1, 'first');
lin_sig_ming_data = [gams(idx_csmg_lin), GM_s1(idx_csmg_lin),...
                     PM_s1(idx_csmg_lin), Sens_gain1(idx_csmg_lin)];
idx_csmg_mpc = find(gams == 0.5, 1, 'first');
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
idx_czmg_lin = find(gams == 2.9, 1, 'first');
lin_zet_ming_data = [gams(idx_czmg_lin), GM_s2(idx_czmg_lin),...
                     PM_s2(idx_czmg_lin), Sens_gain2(idx_czmg_lin)];
idx_czmg_mpc = find(gams == 0.2, 1, 'first');
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
height = 3;
f1 = mkfig(1, width, height); clf
% f1 = figure; clf
freqs = logspace(log10(1), log10(12500), 200);
hsens_cs_rg = frf_bode_mag(Sens_cs_rg, freqs, f1, 'Hz', 'Color', 'k', 'LineStyle', '-');
hsens_cs_rg.DisplayName = 'CS-RG (SLF/MPC)';

hsens_cs_mg_lin = frf_bode_mag(Sens_cs_mg_lin, freqs, f1, 'Hz', 'Color', 'k', 'LineStyle', '-.');
hsens_cs_mg_lin.DisplayName = 'CS-MG (SLF)';
hsens_cs_mg_mpc = frf_bode_mag(Sens_cs_mg_mpc, freqs, f1, 'Hz', 'Color', 'k', 'LineStyle', '--');
hsens_cs_mg_mpcDisplayName = 'CS-MG (MPC)';

hsens_cz_rg = frf_bode_mag(Sens_cz_rg, freqs, f1, 'Hz', 'Color', 'r', 'LineStyle', '-');
hsens_cz_rg.DisplayName = 'CZ-RG (SLF/MPC)';

hsens_cz_mg_lin = frf_bode_mag(Sens_cz_mg_lin, freqs, f1, 'Hz', 'Color', 'r', 'LineStyle', '-.');
hsens_cz_mg_lin.DisplayName = 'CZ-MG (SLF)';
hsens_cz_mg_mpc = frf_bode_mag(Sens_cz_mg_mpc, freqs, f1, 'Hz', 'Color', 'r', 'LineStyle', '--');
hsens_cz_mg_mpc.DisplayName = 'CZ-MG (MPC)';

% response of ppure derivitive.
% frf_bode_mag(zpk([1], [], 1, G.Ts), freqs, f1, 'Hz')

leg_sens = legend([hsens_cs_rg, hsens_cs_mg_lin, hsens_cs_mg_mpc,...
  hsens_cz_rg, hsens_cz_mg_lin, hsens_cz_mg_mpc]);

ax = gca();
ax.XTick = [1          10         100        1000       10000];

tighten_axis(f1, ax)

set(leg_sens, 'FontSize', 7, 'NumColumns', 2, 'Box', 'off');
set(leg_sens, 'Position', [0.5805    0.3410    2.7610    0.4021]);

if saveon
  saveas(f1, fullfile(PATHS.jfig, 'sens_bode.svg'));
end

%%
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

% -------------- Plot Colors & Linestyles ----------------------
GM_colr = 'k';
PM_colr = 'r';
GM_ls1 = '-';
GM_ls2 = '--';

PM_ls1 = '-.';
PM_ls2 = ':';

width = 3.45;
height = 3;
figbase = 10;

f4 = mkfig(4+figbase, width, height); clf
f5 = mkfig(5+figbase, width, height); clf

figure(f4);
ax1 = axes('Position', [0.1100 0.1300 0.7750 0.8150], 'YAxisLocation', 'left');

[hgm1, hpm1] = plot_gmpm(gams, GM_s1, PM_s1, ax1, f4, 'GM_color', GM_colr,...
  'GM_LS', GM_ls1, 'PM_color', PM_colr, 'PM_LS', PM_ls1);
hgm1.DisplayName = 'GM (constant-$\sigma$)';
hpm1.DisplayName = 'PM (constant-$\sigma$)';


figure(f4);
[hgm2, hpm2] = plot_gmpm(gams, GM_s2, PM_s2, ax1, f4, 'GM_color', GM_colr,...
  'GM_LS', GM_ls2, 'PM_color', PM_colr, 'PM_LS', PM_ls2);
hgm2.DisplayName = 'GM (choose-$\zeta$)';
hpm2.DisplayName = 'PM (choose-$\zeta$)';

set(ax1, 'XTick', [0.01, 0.1, 1, 10, 100])
leg1 = legend([hgm1, hgm2, hpm1, hpm2]);
set(leg1, 'Location', 'NorthWest', 'Box', 'off')
%


figure(f5);
ax2 = axes('Position', [0.1100 0.1300 0.7750 0.8150], 'YAxisLocation', 'left');
[h_sens1, h_ts1] = plot_sens_ts(gams, Sens_gain1, TS_s1, ax2, f5,...
  'Ts_color', PM_colr, 'Ts_LS', PM_ls1, 'Sens_color', GM_colr, 'Sens_LS', GM_ls1);
h_sens1.DisplayName = 'Gain of $\mathcal{S}(z)/(z-1)|_{z=1}$ (constant-$\sigma$)';
h_ts1.DisplayName = 'Nominal settle-time (constant-$\sigma$)';

[h_sens2, h_ts2] = plot_sens_ts(gams, Sens_gain2, TS_s2, gca(), f5,...
  'Ts_color', PM_colr, 'Ts_LS', PM_ls2, 'Sens_color', GM_colr, 'Sens_LS', GM_ls2);
h_sens2.DisplayName = 'Gain of $\mathcal{S}(z)/(z-1)|_{z=1}$ (choose-$\zeta$)';
h_ts2.DisplayName = 'Nominal settle-time (choose-$\zeta$)';

leg2 = legend([h_sens1, h_sens2, h_ts1, h_ts2]);
set(leg2, 'Position', [0.1152 0.4316 0.7184 0.2091], 'Box', 'off')
set(ax2, 'XTick', [0.01, 0.1, 1, 10, 100])
yyaxis(ax2, 'left')
ylim([32, 52])
%
figure(f5)
% idx_csro, idx_czro, idx_csmg_lin, idx_csmg_mpc, idx_czmg_lin, idx_czmg_mpc
plot(gams(idx_csro), Sens_gain1(idx_csro), 'xk')
plot(gams(idx_czro), Sens_gain2(idx_czro), 'xk')

plot(gams(idx_csmg_lin), Sens_gain1(idx_csmg_lin), 'ok')
plot(gams(idx_csmg_mpc), Sens_gain1(idx_csmg_mpc), 'sk')

plot(gams(idx_czmg_lin), Sens_gain2(idx_czmg_lin), 'ok')
plot(gams(idx_czmg_mpc), Sens_gain2(idx_czmg_mpc), 'sk')

legend([h_sens1, h_sens2, h_ts1, h_ts2]);

figure(f4);
f4.CurrentAxes = ax1;
yyaxis(ax1, 'left')
plot(gams(idx_csro), GM_s1(idx_csro), 'xk')
yyaxis(ax1, 'right')
plot(gams(idx_csro), PM_s1(idx_csro), 'xk')

yyaxis(ax1, 'left')
plot(gams(idx_csmg_lin), GM_s1(idx_csmg_lin), 'ok')
yyaxis(ax1, 'right')
plot(gams(idx_csmg_lin), PM_s1(idx_csmg_lin), 'ok')

yyaxis(ax1, 'left')
plot(gams(idx_csmg_mpc), GM_s1(idx_csmg_mpc), 'sk')
yyaxis(ax1, 'right')
plot(gams(idx_csmg_mpc), PM_s1(idx_csmg_mpc), 'sk')

% choose zeta
yyaxis(ax1, 'left')
plot(gams(idx_czro), GM_s2(idx_czro), 'xk')
yyaxis(ax1, 'right')
plot(gams(idx_czro), PM_s2(idx_czro), 'xk')

yyaxis(ax1, 'left')
plot(gams(idx_czmg_lin), GM_s2(idx_czmg_lin), 'ok')
yyaxis(ax1, 'right')
plot(gams(idx_czmg_lin), PM_s2(idx_czmg_lin), 'ok')

yyaxis(ax1, 'left')
 plot(gams(idx_czmg_mpc), GM_s2(idx_czmg_mpc), 'sk')
yyaxis(ax1, 'right')
plot(gams(idx_czmg_mpc), PM_s2(idx_czmg_mpc), 'sk')

leg1 = legend([hgm1, hgm2, hpm1, hpm2]);


if saveon
   saveas(f4, fullfile(PATHS.jfig, pmgm_figfile))
   saveas(f5, fullfile(PATHS.jfig, sens_ts_figfile))
end

%%



% ------------------------------------------------------------------
% ------------------------------------------------------------------
% ------- Function Defs. So nice matlab FINALLY allows this --------

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
  ylabel('Gain of $\mathcal{S}(z)/(z-1)|_{z=1}$ [dB]')
  
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




