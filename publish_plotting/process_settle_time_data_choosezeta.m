% This script processes the output data from a sequence of step
% inputs from different experiments (MPC, linear, PID etc) and
% calculates settle-time for each step. The settle-times for each
% experiment are built up into a latex table and saved to a file. 

clear, clc
TOL = 14/512;
tol_mode = 'abs';
verbose = 0;
saveon = true;

addpath(fullfile(getMatPath(), 'afm_mpc_journal', 'functions'))
% where the different experiments are stored.
root = fullfile(PATHS.exp, 'step-exps', 'many_steps_data_rand_18-Jul-2018_01');
% Reference Data with the sequence of references:
load(fullfile(PATHS.exp, 'step-exps', 'many_steps_data_rand_ymax7.mat'))
% reft_pi = load(fullfile(root, 'many_steps_rand.mat'))

whos

Fig = figure(1000); clf
step_ref.plot(Fig);
step_ref.plot_settle_boundary(Fig, TOL, tol_mode);

% ----------------------------------------------------------------
% --------- Constant sigma, rob-opt -------------------
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

names_const_sig_rob_opt = {'LS-CZMG','MPCZ-CZMG',...
                           'LS-CZRO','MPCZ-CZRO',...
                           'LE-CZMG', 'MPCE-CZMG',...
                           'LE-CZRO', 'MPCE-CZRO'};
                  
clrs = {'b', 'r', 'g', 'm', 'k',  [0    0.4470    0.7410],...
    [0.8500    0.3250    0.0980],...
    [0.9290    0.6940    0.1250],...
    [0.4940    0.1840    0.5560],...
    [0.4660    0.6740    0.1880],...
    [0.3010    0.7450    0.9330],...
    [0.6350    0.0780    0.1840]};
  
  
hands = gobjects(1,length(names_const_sig_rob_opt));
TS_s_cell = {};

for k=1:length(names_const_sig_rob_opt)
  dat = load(fullfile(root, files_choose_zet{k}));
  exp_name_str = fields(dat);
  exp_name_str = exp_name_str{1};
  dat = dat.(exp_name_str);
  
%   get_many_steps_ts(dat.y, dat.step_ref.step_amps, dat.step_idx, TOL, verbose, 1, tol_mode);
  ts_s_k = dat.settle_time(TOL, tol_mode, verbose);
  TS_dat_tmp.ts_s = ts_s_k;
  if k==2
%   keyboard
  end

  TS_dat_tmp.name = names_const_sig_rob_opt{k};
  TS_dat_cell{k} = TS_dat_tmp;


%   dat.ploty(Fig, 'Color', clrs{k}, 'LineStyle', '-')
  hands(k) = plot(dat.y.Time, dat.y.Data, 'Color', clrs{k}, 'LineStyle', '-');
  hands(k).DisplayName = TS_dat_tmp.name;
  hold on
end



% ------------------  PI-control --------------------------
% load(fullfile(root,'many_steps_pi.mat'))
% TS_pi = get_many_steps_ts(y_exp, ref_s_pi, step_idx_pi, TOL, verbose, 1, tol_mode);
% TS_dat_tmp.ts_s = TS_pi;
% TS_dat_tmp.name = 'PI';
% TS_dat_cell{end+1} = TS_dat_tmp;

% h5 = plot(y_exp.Time, y_exp.Data, '-m');
% h5.DisplayName = 'PI';

%
clc
% leg = legend([h1, h2, h3,h4, h5]);
leg = legend(hands);


set(leg, 'Location', 'NorthEast')

% ------------------------------------------------------
% ------------ Build the LaTex table -------------------
ts_vec = ManyStepExps.ts_vec_from_dir(root, TOL, tol_mode);
S = ManyStepExps.TS_dat2tex(TS_dat_cell, step_ref, 'do_color', true, 'ts_vec', ts_vec);
fprintf('%s', S);

if saveon
    fid = fopen(fullfile(PATHS.MPCJ_root, 'latex', 'manystepsdata_choosesig.tex'), 'w+');
    fprintf(fid, '%s', S);
    fclose(fid);
end