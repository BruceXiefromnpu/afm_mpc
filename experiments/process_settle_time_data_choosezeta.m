% This script processes the output data from a sequence of step
% inputs from different experiments (MPC, linear, PID etc) and
% calculates settle-time for each step. The settle-times for each
% experiment are built up into a latex table and saved to a file. 

clear, clc

% where the different experiments are stored.
root = fullfile(PATHS.exp, 'step-exps', 'many_steps_data_rand_17-Jul-2018_01');


% Data with the sequence of references:
% load(fullfile(root, 'many_steps_short.mat'))
% reft_pi = load(fullfile(root, 'many_steps.mat'))

load(fullfile(PATHS.exp, 'step-exps', 'many_steps_data_rand_ymax7.mat'))
% reft_pi = load(fullfile(root, 'many_steps_rand.mat'))

saveon = true;

whos
% L = 800;
TOL = 14/512;
tol_mode = 'abs';
verbose = 0;
% ref_s = ref_traj_params.ref_s;
% ref_s_pi = reft_pi.ref_traj_params.ref_s;
% if ref_s ~=ref_s_pi
%   error('need refs the same')
% end

% step_idx = ref_traj_params.impulse_idx;
% step_idx_pi = reft_pi.ref_traj_params.impulse_idx;

Fig = figure(1000); clf
step_ref.plot(Fig)
step_ref.plot_settle_boundary(Fig, TOL, tol_mode);
%%

% load(fullfile(root, 'many_steps_hyst_withsat.mat'))
% TS_hystsat = get_many_steps_ts(y_exp, ref_s, step_idx, TOL, 1);

% h1 = plot(y_exp.Time, y_exp.Data, '-r');
% h1.DisplayName = 'linfp, inv hyst w/ sat ';
%


% many_steps_ymax7_lin_EXP_choose-zet-min-gam_07-17-2018.mat    
% many_steps_ymax7_mpc_EXP_choose-zet-min-gam_07-17-2018.mat
% many_steps_ymax7_linfxp_sim_choose-zet-min-gam_07-17-2018.mat
% many_steps_ymax7_mpcfxp_sim_choose-zet-min-gam_07-17-2018.mat

% many_steps_ymax7_linfxp_sim_const-sig-min-gam_07-17-2018.mat
% many_steps_ymax7_mpcfxp_sim_const-sig-min-gam_07-17-2018.mat
% many_steps_ymax7_lin_EXP_const-sig-rob-opt_07-17-2018.mat
% many_steps_ymax7_mpc_EXP_const-sig-rob-opt_07-17-2018.mat
% 
% many_steps_ymax7_linfxp_sim_const-sig-rob-opt_07-17-2018.mat
% many_steps_ymax7_mpcfxp_sim_const-sig-rob-opt_07-17-2018.mat
% many_steps_ymax7_lin_EXP_const-sig-min-gam_07-17-2018.mat
% many_steps_ymax7_mpc_EXP_const-sig-min-gam_07-17-2018.mat
% 


clc
TS_s_cell = {};
% step_idx = step_idx(1:end-2);
% ref_s = ref_s(1:end-2);
% ----------------------------------------------------------------
% --------- Constant sigma, rob-opt -------------------
files_const_sig_rob_opt = {
'many_steps_ymax7_linfxp_sim_choose-zet-rob-opt_07-17-2018.mat',...
'many_steps_ymax7_mpcfxp_sim_choose-zet-rob-opt_07-17-2018.mat',...
'many_steps_ymax7_lin_EXP_choose-zet-rob-opt_07-17-2018.mat',...
'many_steps_ymax7_mpc_EXP_choose-zet-rob-opt_07-17-2018.mat',...
'many_steps_ymax7_mpcfxp_sim_choose-zet-min-gam_07-17-2018.mat',...
'many_steps_ymax7_linfxp_sim_choose-zet-min-gam_07-17-2018.mat',...
'many_steps_ymax7_lin_EXP_choose-zet-min-gam_07-17-2018.mat',...
'many_steps_ymax7_mpc_EXP_choose-zet-min-gam_07-17-2018.mat',...
};


names_const_sig_rob_opt = {'LS (CZMG)','MPCS (CZMG)',...
                   'LE (CZRO)', 'MPCE (CZRO)',...
                   'LS (CZMG)','MPCS (CZMG)',...
                    'LE (CZMG)', 'MPCE (CZMG)'};



clrs = {'b', 'r', 'g', 'm', 'k',  [0    0.4470    0.7410],...
    [0.8500    0.3250    0.0980],...
    [0.9290    0.6940    0.1250],...
    [0.4940    0.1840    0.5560],...
    [0.4660    0.6740    0.1880],...
    [0.3010    0.7450    0.9330],...
    [0.6350    0.0780    0.1840]};
  
  
hands = gobjects(1,length(names_const_sig_rob_opt));
for k=1:length(names_const_sig_rob_opt)
  dat = load(fullfile(root, files_const_sig_rob_opt{k}));
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

% -- First, we programmatically construct \tabular{ccc},
%    since the 'ccc' depends on how many columns we need.
%    the number of columns in the table
c_fmt = repmat('c', 1, 3+length(TS_dat_cell)); 
S = sprintf('\\begin{tabular}{%s}\n', c_fmt);

% -- Form the table header:
str_ref_cols = sprintf('&ref & delta');
str_dat_cols = '';
for k=1:length(TS_dat_cell)
    str_dat_cols = sprintf(' %s & %s', str_dat_cols, TS_dat_cell{k}.name);
end
S = sprintf('%s%s%s\\\\\n\\toprule\n', S, str_ref_cols, str_dat_cols)
%
% -- Build up the body of the table. Outer loop is for each row.
% Inner loop is for each experiment (columns).

for k = 2:length(step_ref.step_amps)
  
  delta_ref = step_ref.step_amps(k) - step_ref.step_amps(k-1);
  str_ref_cols = sprintf('&%.2f & %.2f', step_ref.step_amps(k), delta_ref);
  str_dat_cols = '';
  for j = 1:length(TS_dat_cell)
      str_dat_cols = sprintf('%s &%.2f', str_dat_cols, ...
                             1000*TS_dat_cell{j}.ts_s(k-1));
      % keyboard
  end
  s_row = sprintf('%s%s\\\\ \n', str_ref_cols, str_dat_cols);
  S = sprintf('%s%s', S, s_row);
end
% -- Build the footer. This is where the totals go.
str_ref_cols = sprintf('total & -- & --');
str_dat_cols = '';
for j = 1:length(TS_dat_cell)
    str_dat_cols = sprintf('%s &%.2f', str_dat_cols, ...
                           1000*sum(TS_dat_cell{j}.ts_s) );
end

s_row = sprintf('%s%s\\\\ \n', str_ref_cols, str_dat_cols);

S = sprintf('%s\\midrule\n %s', S, s_row);
               

S = sprintf('%s\\end{tabular}\n', S)


if saveon
    fid = fopen(fullfile(PATHS.MPCJ_root, 'latex', 'manystepsdata_choosesig.tex'), 'w+');
    fprintf(fid, '%s', S);
    fclose(fid);
end