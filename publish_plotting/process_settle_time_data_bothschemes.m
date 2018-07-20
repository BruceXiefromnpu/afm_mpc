% This script processes the output data from a sequence of step
% inputs from different experiments (MPC, linear, PID etc) and
% calculates settle-time for each step. The settle-times for each
% experiment are built up into a latex table and saved to a file. 

clear, clc
saveon = false;
TOL = 14/512;
tol_mode = 'abs';
verbose = 0;

addpath(fullfile(getMatPath(), 'afm_mpc_journal', 'functions'))
% where the different experiments are stored.
root = fullfile(PATHS.exp, 'step-exps', 'many_steps_data_rand_17-Jul-2018_01');

% Reference Data 
load(fullfile(PATHS.exp, 'step-exps', 'many_steps_data_rand_ymax7.mat'))
% reft_pi = load(fullfile(root, 'many_steps_rand.mat'))

whos

Fig = figure(1000); clf
step_ref.plot(Fig);
step_ref.plot_settle_boundary(Fig, TOL, tol_mode);



% ----------------------------------------------------------------
% --------- Constant sigma, rob-opt -------------------
files_const_sig = {
'many_steps_ymax7_linfxp_sim_const-sig-min-gam_07-17-2018.mat',...
'many_steps_ymax7_mpcfxp_sim_const-sig-min-gam_07-17-2018.mat',...
'many_steps_ymax7_linfxp_sim_const-sig-rob-opt_07-17-2018.mat',...  
'many_steps_ymax7_mpcfxp_sim_const-sig-rob-opt_07-17-2018.mat',...
'many_steps_ymax7_lin_EXP_const-sig-min-gam_07-17-2018.mat',...
'many_steps_ymax7_mpc_EXP_const-sig-min-gam_07-17-2018.mat',...
'many_steps_ymax7_lin_EXP_const-sig-rob-opt_07-17-2018.mat',...
'many_steps_ymax7_mpc_EXP_const-sig-rob-opt_07-17-2018.mat'};

names_const_sig_rob_opt = {'LS-CSMG','MPCS-CSMG',...
                           'LS-CSRO','MPCS-CSRO',...
                           'LE-CSMG', 'MPCE-CSMG',...
                            'LE-CSRO', 'MPCE-CSRO'};
% 
% clrs = {[0    0.4470    0.7410],...
%     [0.8500    0.3250    0.0980],...
%     [0.9290    0.6940    0.1250],...
%       [0.4940    0.1840    0.5560],...};
clrs = {'b', 'r', 'g', 'k', 'b', 'r', 'g', 'k'}    ;
line_styles = {'-', '--', '-', '--','-', '--','-', '--','-', '--','-', '--','-', '--'};
    
  


% hands = gobjects(1,length(names_const_sig_rob_opt));
hands = gobjects(1,4);
TS_s_cell = {};

step_exps_CS_cell = cell(1, length(names_const_sig_rob_opt));
for k=1:length(names_const_sig_rob_opt)
  dat = load(fullfile(root, files_const_sig{k}));
  exp_name_str = fields(dat);
  exp_name_str = exp_name_str{1};
  dat = dat.(exp_name_str);
  
  dat.name = names_const_sig_rob_opt{k};
  dat.Color = clrs{k};
  dat.LineStyle = line_styles{k};
  step_exps_CS_cell{k} = dat;

end

% leg = legend(hands);
% set(leg, 'Location', 'NorthEast')

step_exps_CS = ManyStepExps(TOL, tol_mode, step_ref, step_exps_CS_cell{:});

% ------------ Build the LaTex table -------------------
ts_vec = ManyStepExps.ts_vec_from_dir(root, TOL, tol_mode);
S = step_exps_CS.TS_dat2tex('do_color', true, 'ts_vec', ts_vec);
fprintf('%s', S); % just display it.
if saveon
  ManyStepExps.write_tex_data(S, fullfile(PATHS.MPCJ_root, 'latex', 'manystepsdata.tex'));
end

%% Now, we can plot

width = 3.4;
height = 3.5;

Fig = mkfig(10, width, height); 
%
clf, clc
% subplot(2,2,[1,2])
ax1 = axes('Position', [0.1300 0.75 0.7750 0.25], 'Units', 'normalized');
ax2 = axes('Position', [0.1300 0.1100 0.36 0.55], 'Units', 'normalized'); 
ax3 = axes('Position', [0.5703 0.1100 .36 0.55], 'Units', 'normalized');
% xlabel(



step_ref.plot(ax1);
step_ref.plot_settle_boundary(ax1, TOL, tol_mode);
step_exps_CS.ploty_selected([5:8], ax1);
ax1.YLabel.String = 'y [v]';
ax1.XLabel.String = 't [ms]';
set(ax1, 'YLim', [-7.5, 7.5]);

% "best"
annotation('ellipse', [0.18,0.78, 0.04, 0.05], 'Color', 'r');
annotation('arrow', [0.2, 0.2], [0.785, 0.6]); %, 'Color', 'k');
% "worst"
annotation('ellipse', [0.55, 0.93, 0.04, 0.05], 'Color', 'r');
annotation('arrow', [0.57, 0.59], [0.95, 0.6]); %, 'Color', 'k');

% Even though we want to zoom in, its maybe easiest to stick with our framework
% and plot the whole thing, then adjust xlim and ylim.

step_ref.plot(ax2);
step_ref.plot_settle_boundary(ax2, TOL, tol_mode);
step_exps_CS.ploty_selected([5:8], ax2);
set(ax2, 'XLim', [0.20, 0.244]); % ;
set(ax2, 'YLim', [-4.3, -3.8]);
ax2.XLabel.String = 't [ms]';
ax2.YLabel.String = 'y [v]';

step_ref.plot(ax3);
step_ref.plot_settle_boundary(ax3, TOL, tol_mode);
hands = step_exps_CS.ploty_selected([5:8], ax3);

set(ax3, 'XLim', [1.403, 1.429]);
ax3.XLabel.String = 't [ms]';
set(ax3, 'YLim', [5, 5.6]);

leg = legend(hands);
set(leg, 'FontSize', 8, 'Box', 'off',...
  'Position', [0.650, 0.5040, 0.3113, 0.1307])

%%
if saveon
  saveas(Fig, fullfile(PATHS.jfig, 'step_exps_const_sig_y.svg'))
end
%%

%% Now, the choos-zeta scenario
% ----------------------------------------------------------------
% --------- Choose, rob-opt -----------------------------
clear dat; clear step_exps_cell; clear step_exps;
root = fullfile(PATHS.exp, 'step-exps', 'many_steps_data_rand_18-Jul-2018_01');

files_choose_zet = {
'many_steps_ymax7_linfxp_sim_choose-zet-min-gam_07-18-2018.mat',...
'many_steps_ymax7_mpcfxp_sim_choose-zet-min-gam_07-18-2018.mat',...
'many_steps_ymax7_linfxp_sim_choose-zet-rob-opt_07-18-2018.mat',...
'many_steps_ymax7_mpcfxp_sim_choose-zet-rob-opt_07-18-2018.mat',...
'many_steps_ymax7_lin_EXP_choose-zet-min-gam_07-18-2018.mat',...
'many_steps_ymax7_mpc_EXP_choose-zet-min-gam_07-18-2018.mat',...
'many_steps_ymax7_lin_EXP_choose-zet-rob-opt_07-18-2018.mat',...
'many_steps_ymax7_mpc_EXP_choose-zet-rob-opt_07-18-2018.mat',...
};

names_choose_zet = {'LS-CZMG','MPCZ-CZMG',...
                           'LS-CZRO','MPCZ-CZRO',...
                           'LE-CZMG', 'MPCE-CZMG',...
                           'LE-CZRO', 'MPCE-CZRO'};

clrs = {'b', 'r', 'g', 'k', 'b', 'r', 'g', 'k'}    ;
line_styles = {'-', '--', '-', '--','-', '--','-', '--','-', '--','-', '--','-', '--'};

step_exps_CZ_cell = cell(1, length(files_choose_zet));
for k=1:length(names_choose_zet)
  dat = load(fullfile(root, files_choose_zet{k}));
  exp_name_str = fields(dat);
  exp_name_str = exp_name_str{1};
  dat = dat.(exp_name_str);
  
  dat.name = names_choose_zet{k};
  dat.Color = clrs{k};
  dat.LineStyle = line_styles{k};
  step_exps_CZ_cell{k} = dat;

end

% leg = legend(hands);
% set(leg, 'Location', 'NorthEast')

step_exps_CZ = ManyStepExps(TOL, tol_mode, step_ref, step_exps_CZ_cell{:});

% ------------ Build the LaTex table -------------------
ts_vec = ManyStepExps.ts_vec_from_dir(root, TOL, tol_mode);
S = step_exps_CZ.TS_dat2tex('do_color', true, 'ts_vec', ts_vec);
fprintf('%s', S); % just display it.
if saveon
  ManyStepExps.write_tex_data(S, fullfile(PATHS.MPCJ_root, 'latex', 'manystepsdata_choosezeta.tex'));
end

%% Now, we can plot

width = 3.4;
height = 3.5;

Fig = mkfig(11, width, height); 
%
clf, clc
% subplot(2,2,[1,2])
ax1 = axes('Position', [0.1300 0.75 0.7750 0.25], 'Units', 'normalized');
ax2 = axes('Position', [0.1300 0.1100 0.36 0.55], 'Units', 'normalized'); 
ax3 = axes('Position', [0.5703 0.1100 .36 0.55], 'Units', 'normalized');
% xlabel(



step_ref.plot(ax1);
step_ref.plot_settle_boundary(ax1, TOL, tol_mode);
step_exps_CZ.ploty_selected([5:8], ax1);
ax1.YLabel.String = 'y [v]';
ax1.XLabel.String = 't [ms]';
set(ax1, 'YLim', [-7.5, 7.5]);

% "best"
annotation('ellipse', [0.18,0.78, 0.04, 0.05], 'Color', 'r');
annotation('arrow', [0.2, 0.2], [0.785, 0.6]); %, 'Color', 'k');
% "worst"
annotation('ellipse', [0.55, 0.93, 0.04, 0.05], 'Color', 'r');
annotation('arrow', [0.57, 0.59], [0.95, 0.6]); %, 'Color', 'k');

% Even though we want to zoom in, its maybe easiest to stick with our framework
% and plot the whole thing, then adjust xlim and ylim.

step_ref.plot(ax2);
step_ref.plot_settle_boundary(ax2, TOL, tol_mode);
step_exps_CZ.ploty_selected([5:8], ax2);
set(ax2, 'XLim', [0.20, 0.244]); % ;
set(ax2, 'YLim', [-4.3, -3.8]);
ax2.XLabel.String = 't [ms]';
ax2.YLabel.String = 'y [v]';

step_ref.plot(ax3);
step_ref.plot_settle_boundary(ax3, TOL, tol_mode);
hands = step_exps_CZ.ploty_selected([5:8], ax3);
%%
set(ax3, 'XLim', [2.1, 2.15]);
ax3.XLabel.String = 't [ms]';
set(ax3, 'YLim', [6.8, 7.2]);

leg = legend(hands);
set(leg, 'FontSize', 8, 'Box', 'off',...
  'Position', [0.650, 0.5040, 0.3113, 0.1307])

%%
if saveon
  saveas(Fig, fullfile(PATHS.jfig, 'step_exps_choose_zet_y.svg'))
end




