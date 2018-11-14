% This script processes the output data from a sequence of step
% inputs from different experiments (MPC, linear, PID etc) and
% calculates settle-time for each step. The settle-times for each
% experiment are built up into a latex table and saved to a file. 

% The ostensible difference between this scripts (**_v3.m) and (**_v2.m is
% is to plot the min-gamma on the same figure and rob-gamma on a different
% figure, rather than CS on one figure and CZ on the other figure. 

clear, clc
saveon = true;
TOL = 14/512;
tol_mode = 'abs';
verbose = 0;

addpath(fullfile(getMatPath(), 'afm_mpc_journal', 'functions'))
% where the different experiments are stored.

% Reference Data 
load(fullfile(PATHS.exp, 'step-exps', 'many_steps_data_rand_ymax7_n6p5.mat'))
% reft_pi = load(fullfile(root, 'many_steps_rand.mat'))

whos

Fig = figure(1000); clf
step_ref.yscaling = 5;
step_ref.plot(Fig);
step_ref.plot_settle_boundary(Fig, TOL, tol_mode);
ax0 = gca();
%%
clc


files_min_gam = {
'many_steps_ymax7_linfxp_sim_const-sig_09-21-2018_gam_7.50.mat',...
'many_steps_ymax7_mpcfxp_sim_const-sig_09-21-2018_gam_0.00.mat',...
'many_steps_ymax7_linfxp_sim_choose-zet_09-21-2018_gam_3.50.mat',...
'many_steps_ymax7_mpcfxp_sim_choose-zet_09-21-2018_gam_0.00.mat',...
'many_steps_ymax7_lin_EXP_const-sig_09-24-2018_gam_7.50_obsID_3.mat',...
'many_steps_ymax7_mpc_EXP_const-sig_09-24-2018_gam_0.00_obsID_3.mat',...
'many_steps_ymax7_lin_EXP_choose-zet_09-24-2018_gam_3.50_obsID_3.mat',...
'many_steps_ymax7_mpc_EXP_choose-zet_09-24-2018_gam_0.00_obsID_3.mat',...
};

gam_str_cs = 'gam_100.00'
gam_str_cz = 'gam_25.00'
files_rob_gam = {
  ['many_steps_ymax7_linfxp_sim_const-sig_09-21-2018_', gam_str_cs, '.mat'],...
['many_steps_ymax7_mpcfxp_sim_const-sig_09-21-2018_', gam_str_cs, '.mat'],...
['many_steps_ymax7_linfxp_sim_choose-zet_09-21-2018_', gam_str_cz, '.mat'],...
['many_steps_ymax7_mpcfxp_sim_choose-zet_09-21-2018_', gam_str_cz, '.mat'],...
['many_steps_ymax7_lin_EXP_const-sig_09-24-2018_', gam_str_cs, '_obsID_3.mat'],...
['many_steps_ymax7_mpc_EXP_const-sig_09-24-2018_', gam_str_cs, '_obsID_3.mat'],...
['many_steps_ymax7_lin_EXP_choose-zet_09-24-2018_', gam_str_cz, '_obsID_3.mat'],...
['many_steps_ymax7_mpc_EXP_choose-zet_09-24-2018_',gam_str_cz, '_obsID_3.mat'],...
};

% Indeces 1-4 are simultion.
cs_idx = [5,6];
cz_idx = [7,8];

names_min_gam = {
  'SLF-CS ($\gamma=7.5$)','MPC-CS ($\gamma=10^{-3}$',...
  'SLF-CZ-($\gamma=3.5$)', 'MPC-CZ ($\gamma=10^{-5})$',...
  'SLF-CS ($\gamma=7.5$)','MPC-CS ($\gamma=10^{-3}$)',...
  'SLF-CZ-($\gamma=3.5$)', 'MPC-CZ ($\gamma=10^{-5}$)'};


names_rob_gam = {
  'SLF-CS ($\gamma=100)','MPC-CS ($\gamma=100$)',...
  'SLF-CZ ($\gamma=25)', 'MPC-CZ ($\gamma=25$)',...
  'SLF-CS ($\gamma=100$)', 'MPC-CS ($\gamma=100$)',...
  'SLF-CZ ($\gamma=25$)', 'MPC-CZ ($\gamma=25$)'};

data_root = fullfile(PATHS.exp, 'step-exps', 'many_steps_sweep_gamma_21-Sep-2018_01');
clrs = {'b', 'r', 'g', 'k', 'b', 'r', 'g', 'k'}    ;
line_styles = {'-', '--', '-', '--','-', '--','-', '--','-', '--','-', '--','-', '--'};

% ----------------------------------------------------------------
% --------- Load min-gamma data --------- -------------------
step_exps_MG_cell = cell(1, length(names_min_gam));
for k=1:length(names_min_gam)
  dat = load(fullfile(data_root, files_min_gam{k}));
  exp_name_str = fields(dat);
  exp_name_str = exp_name_str{1};
  dat = dat.(exp_name_str);
  
  dat.Ipow = dat.Ipow * 1000/15.15;
  dat.yscaling = 5;
  dat.yunits = '[$\mu$m]';

  dat.name = names_min_gam{k};
  dat.Color = clrs{k};
  dat.LineStyle = line_styles{k};
  step_exps_MG_cell{k} = dat;
end

step_exps_MG = ManyStepExps(TOL, tol_mode, step_ref, step_exps_MG_cell{:});



% -------------------------------------------------------------------------
% ------------- Load robust-gamma data -------------------------------------

% clrs = {'b', 'r', 'g', 'k', 'g', 'k', 'b', 'r'};
clrs = {'b', 'r', 'g', 'k', 'b', 'r', 'g', 'k'};
line_styles = {'-', '--', '-', '--','-', '--','-', '--','-', '--','-', '--','-', '--'};

step_exps_RG_cell = cell(1, length(files_rob_gam));
for k=1:length(names_rob_gam)
  dat = load(fullfile(data_root, files_rob_gam{k}));
  
  exp_name_str = fields(dat);
  exp_name_str = exp_name_str{1};
  dat = dat.(exp_name_str);
  dat.yscaling = 5;
  dat.yunits = '[$\mu$m]';
  dat.Ipow = dat.Ipow * 1000/15.15;
  dat.name = names_rob_gam{k};
  dat.Color = clrs{k};
  dat.LineStyle = line_styles{k};
  step_exps_RG_cell{k} = dat;

end

step_exps_RG = ManyStepExps(TOL, tol_mode, step_ref, step_exps_RG_cell{:});

step_exps_MG.ploty_selected(cs_idx, ax0);
step_exps_RG.ploty_selected(cs_idx, ax0);

%%
ts_master_vec = unique([step_exps_RG.TS_mat(:); step_exps_MG.TS_mat(:)]);

% -------------------------------------------------------------------------
rgb1 = [0.230, 0.299, 0.754];
rgb2 = [0.706, 0.016, 0.150];
s_ = linspace(0,1, length(step_ref.step_diff_amps));
color_map = diverging_map(s_, rgb1, rgb2);
% ------------ Constant-Sigma LaTex table -------------------
% ts_master_vec = ManyStepExps.ts_vec_from_dir(data_root, TOL, tol_mode);

S = step_exps_MG.TS_dat2tex('do_color', true, 'ts_vec', ts_master_vec, 'colormap', color_map);
fprintf('%s', S); % just display it.

if saveon
  ManyStepExps.write_tex_data(S, fullfile(PATHS.MPCJ_root, 'latex', 'manystepsdata_mingam.tex'));
end


% ------------ Choose-Zeta LaTex table -------------------
S = step_exps_RG.TS_dat2tex('do_color', true, 'ts_vec', ts_master_vec, 'colormap', color_map);
fprintf('%s', S); % just display it.
if saveon
  ManyStepExps.write_tex_data(S, fullfile(PATHS.MPCJ_root, 'latex', 'manystepsdata_choosezeta.tex'));
end

% create a figure which is only a colormap legend
fig100 = mkfig(100, 7, 0.75); clf
ax = gca();
set(ax, 'Visible', 'off');
colormap(ax, color_map);
cb = colorbar(ax);
set(cb, 'Position', [0.03, 0.5, 0.94, 0.25], 'Orientation', 'horizontal',...
  'Units', 'normalized', 'AxisLocation', 'in', 'FontSize', 9)
cb.Label.String = 'settle-time [ms]';
ts_min = min(ts_master_vec)*1000;
ts_max = max(ts_master_vec)*1000;
caxis([ts_min, ts_max+.001]); % +.001 to get 93 to display
set(cb, 'Ticks', [3, 25, 50, 75, 93])
saveas(fig100, fullfile(PATHS.jfig(), 'ts_colorbar.svg'))


%% -------------------------------------------------------------------------
% ------------ Minimum-gamma Plot steps, experimental-only zoom in--------
figure(1000); clf
ax0 = gca();
step_ref.plot(ax0);
step_ref.plot_settle_boundary(ax0, TOL, tol_mode);
step_exps_MG.ploty_selected([cs_idx, cz_idx], ax0);

%%
clc
width = 3.4;
height = 3.75;

Fig = mkfig(10, width, height); clf

ax1 = axes('Units', 'inches', 'Position', [0.35, 0.3850+1.25, 1.35, 1.9250]);
ax2 = axes('Units', 'inches', 'Position', [1.99, 0.3850+1.25, 1.35, 1.9250]);

ax3 = axes('Units', 'inches', 'Position', [0.35, 0.385, 1.35, 1.0], 'Box', 'on');
ax4 = axes('Units', 'inches', 'Position', [1.99, 0.385, 1.35, 1.0], 'Box', 'on');

% Even though we want to zoom in, its maybe easiest to stick with our framework
% and plot the whole thing, then adjust xlim and ylim.
xlim_lft = [0.50, 0.53];
xlim_rt = [2.30, 2.33];
ylim_lft = [7, 8.5];
ylim_rt = [-1.5, 2.1];

step_ref.plot(ax1);
step_ref.plot_settle_boundary(ax1, TOL, tol_mode);
step_exps_MG.ploty_selected([cs_idx, cz_idx], ax1);

set(ax1, 'XLim', xlim_lft); % ;
set(ax1, 'YLim', ylim_lft);
ax1.YLabel.String = 'y [$\mu$m]';
title(ax1, 'step 5')

step_ref.plot(ax2);
step_ref.plot_settle_boundary(ax2, TOL, tol_mode);
hands = step_exps_MG.ploty_selected([cs_idx, cz_idx], ax2);

set(ax2, 'XLim', xlim_rt);
set(ax2, 'YLim', ylim_rt);
title(ax2, 'step 24')

leg = legend(hands);
set(leg, 'FontSize', 7, 'Box', 'off',...
  'Position',  [0.5811 0.8110 0.4121 0.1397]);


hold(ax3, 'on')
step_exps_MG.plotdu_selected([cs_idx, cz_idx], ax3);

set(ax3, 'XLim', xlim_lft-[0, 0.025]);
grid(ax3, 'on')
ax3.YLabel.String='$\Delta u_k$';

hold(ax4, 'on')
step_exps_MG.plotdu_selected([cs_idx, cz_idx], ax4);

set(ax4, 'XLim', xlim_rt-[0, 0.025]);
set(ax4, 'XLim', xlim_rt-[0, 0.025]);
grid(ax4, 'on')
ax3.YLabel.String = '$\Delta u_k$';

xlabel(ax3, 'time [s]')
xlabel(ax4, 'time [s]')
ylim(ax4, [-0.2, 0.2])
set(ax1, 'XTick', [0.5, 0.51, 0.52])
set(ax2, 'XTick', [2.3, 2.31, 2.32])
set(ax3, 'XTick', [0.5, 0.5025, 0.505])
set(ax4, 'YTick', [-0.1, 0, 0.1, 0.2]);

if saveon
  saveas(Fig, fullfile(PATHS.jfig, 'step_exps_min_gam.svg'))
end
%% Now, the choos-zeta scenario
% ----------------------------------------------------------------
Fig = mkfig(1001, width, 1.5); clf
ax00 = gca();
step_ref.plot(ax00);
xlabel(ax00, 'time [s]')
ylabel('reference [$\mu$m]')
tighten_axis(Fig, ax00)

% text(0.5, 11, '5')
% text(2.3, 5, '24')
a1 = annotation('ellipse');
set(a1, 'Units', 'inches', 'Position', [0.9245 0.8966 0.1360 0.1750], 'Color', 'r');
a3 = annotation('ellipse');
set(a3, 'Units', 'inches', 'Position', [3.0339 0.7917 0.1640 0.1562], 'Color', 'r')


saveas(Fig, fullfile(PATHS.jfig, 'steps.svg'))
%%
width = 3.4;
height = 3.75;
Fig = mkfig(11, width, height); clf
%
clf, clc
ax1 = axes('Units', 'inches', 'Position', [0.35, 0.3850+1.25, 1.35, 1.9250]);
ax2 = axes('Units', 'inches', 'Position', [1.99, 0.3850+1.25, 1.35, 1.9250]);

ax3 = axes('Units', 'inches', 'Position', [0.35, 0.385, 1.35, 1.0], 'Box', 'on');
ax4 = axes('Units', 'inches', 'Position', [1.99, 0.385, 1.35, 1.0], 'Box', 'on');
% Even though we want to zoom in, its maybe easiest to stick with our framework
% and plot the whole thing, then adjust xlim and ylim.

xlim_lft = [0.50, 0.53];
xlim_rt = [2.30, 2.33];
ylim_lft = [7, 8.5];
ylim_rt = [-2, 2];

step_ref.plot(ax1);
step_ref.plot_settle_boundary(ax1, TOL, tol_mode);
step_exps_RG.ploty_selected([cs_idx, cz_idx], ax1);

set(ax1, 'XLim', [0.50, 0.53]); % ;
set(ax1, 'YLim', ylim_lft);
set(ax1, 'XTick', [0.5, 0.51, 0.52])

ax1.YLabel.String = 'y [$\mu$m]';
title(ax1, 'step 5')
step_ref.plot(ax2);
step_ref.plot_settle_boundary(ax2, TOL, tol_mode);
hands = step_exps_RG.ploty_selected([cs_idx, cz_idx], ax2);

set(ax2, 'XLim', xlim_rt);
set(ax2, 'YLim', ylim_rt);
set(ax2, 'XTick', [2.3, 2.31, 2.32])
title(ax2, 'step 24')
leg = legend(hands);
set(leg, 'FontSize', 7, 'Box', 'off',...
  'Position', [0.5882 0.7677 0.3978 0.1397]);

hold(ax3, 'on')
step_exps_RG.plotdu_selected([cs_idx, cz_idx], ax3);
set(ax3, 'XLim', xlim_lft-[0, 0.025]);
grid(ax3, 'on')
ax3.YLabel.String='$\Delta u_k$';
set(ax3, 'XTick', [0.5, 0.5025, 0.505])

hold(ax4, 'on')
step_exps_RG.plotdu_selected([cs_idx, cz_idx], ax4);
set(ax4, 'XLim', xlim_rt-[0, 0.025]);
% set(ax5, 'XTick', [2.1, 2.125, 2.15])
grid(ax4, 'on')

xlabel(ax3, 'time [s]')
xlabel(ax4, 'time [s]')

ylim(ax4, [-0.2, 0.2])
set(ax4, 'YTick', [-0.1, 0, 0.1, 0.2]);

if saveon
  saveas(Fig, fullfile(PATHS.jfig, 'step_exps_rob_gam.svg'))
end

%%

 clc

data_root = fullfile(PATHS.exp, 'step-exps', 'many_steps_sweep_gamma_21-Sep-2018_01');
files = split(ls(data_root));

Imax_s = [];

for k =1:length(files)
  file = files{k};
  
  if ~contains(file, 'obsID_')
    disp(file)
    continue
  end
  
  if mod(k, 10) == 0
    fprintf('loading file number %d\n', k)
  end
  
  dat = load(fullfile(data_root, file));
  
  if isfield(dat, 'afm_exp_lin')
    field_name = 'afm_exp_lin';
  elseif isfield(dat, 'afm_exp_mpc')
    field_name = 'afm_exp_mpc';
  end
  
  ipow_data = dat.(field_name).Ipow;
  Imax_s(k) = max(ipow_data.Data* 1000/15.15);
%   keyboard
end

fprintf('max of all currents: %.4f [mA]\n', max(Imax_s))
%%


% Now, plot current draw
clc
height = 2.5;
F12 = mkfig(12, width, height);
ax = gca();

hands_Ipow = step_exps_CS.plotIpow_selected(5:8, ax);

ylim([-130, 110])
grid on
ylabel('Current [mA]')
xlabel('time [s]')

leg3 = legend(hands_Ipow);
set(leg3, 'NumColumns', 2, 'Units', 'inches',...
  'Position', [0.8089 0.2993 2.5049 0.3319], 'Box', 'off')
tighten_axis(F12, ax)


F13 = mkfig(13, width, height);
ax = gca();

hands_Ipow = step_exps_RG.plotIpow_selected(5:8, ax);

ylim([-130, 110])
grid on
ylabel('Current [mA]')
xlabel('time [s]')

leg4 = legend(hands_Ipow);
set(leg4, 'NumColumns', 2, 'Units', 'inches',...
  'Position', [0.7956 0.2993 2.5181 0.3319], 'Box', 'off')
tighten_axis(F13, ax)



%%
saveas(F12, fullfile(PATHS.jfig, 'step_exps_const_sig_Ipow.svg'))
saveas(F13, fullfile(PATHS.jfig, 'step_exps_choose_zet_Ipow.svg'))
