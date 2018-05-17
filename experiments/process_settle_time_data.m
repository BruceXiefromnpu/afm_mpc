% This script processes the output data from a sequence of step
% inputs from different experiments (MPC, linear, PID etc) and
% calculates settle-time for each step. The settle-times for each
% experiment are built up into a latex table and saved to a file. 

clear, clc

% where the different experiments are stored.
root = fullfile(PATHS.exp, 'experiments', 'many_steps_data')
% addpath('many_steps_data')

% Data with the sequence of references:
load('many_steps.mat')

saveon = true;

whos
L = 800;
TOL = 0.01;
tol_mode = 'rel';
verbose = 0;
ref_s = ref_traj_params.ref_s;
step_idx = ref_traj_params.impulse_idx;


figure(1000); clf
plot(ref_traj_params.ref_traj.Time, ref_traj_params.ref_traj.Data)
hold on, grid on;

r = ref_traj_params.ref_traj;
for k=2:length(ref_s)
    idx_start = step_idx(k);
    if k == length(ref_s)
      idx_end = length(r.Time);
    else
      idx_end = step_idx(k+1)-1;
    end

    ref = ref_s(k);
    
    t0 = r.Time(idx_start);
    t1 = r.Time(idx_end);
    
    plot([t0, t1], [ref+TOL, ref+TOL], ':k');
    plot([t0, t1], [ref-TOL, ref-TOL], ':k');
end



% load(fullfile(root, 'many_steps_hyst_withsat.mat'))
% TS_hystsat = get_many_steps_ts(y_exp, ref_s, step_idx, TOL, 1);

% h1 = plot(y_exp.Time, y_exp.Data, '-r');
% h1.DisplayName = 'linfp, inv hyst w/ sat ';
%

TS_s_cell = {};

load(fullfile(root,'many_steps_hyst_nosat_R4.mat'))
TS_hyst = get_many_steps_ts(y_exp, ref_s, step_idx, TOL, verbose, 1, tol_mode);
TS_dat_tmp.ts_s = TS_hyst;
TS_dat_tmp.name = 'Lin (fp) w/ $\mathcal{H}^{-1}$';
TS_dat_cell{1} = TS_dat_tmp;

h2 = plot(y_exp.Time, y_exp.Data, '-b');
h2.DisplayName = 'Invert hyst (linfp, no sat operator)';


% ------------- Linear, with no inversion -----------------
load(fullfile(root,'many_steps_noinvert.mat'))
TS_lin_noinv = get_many_steps_ts(y_exp, ref_s, step_idx, TOL, verbose, 1, tol_mode);
TS_dat_tmp.ts_s = TS_lin_noinv;
TS_dat_tmp.name = 'Lin (fp) w/ no-inv';
TS_dat_cell{end+1} = TS_dat_tmp;

h3 = plot(y_exp.Time, y_exp.Data, '-g');
h3.DisplayName = 'linfp, no inversion';

% ------------------  PI-control --------------------------
load(fullfile(root,'many_steps_pi.mat'))
TS_pi = get_many_steps_ts(y_exp, ref_s, step_idx, TOL, verbose, 1, tol_mode);
TS_dat_tmp.ts_s = TS_pi;
TS_dat_tmp.name = 'PI';
TS_dat_cell{end+1} = TS_dat_tmp;

h4 = plot(y_exp.Time, y_exp.Data, '-m');
h4.DisplayName = 'PI';


load(fullfile(root,'many_steps_mpc_invHyst_invDrift.mat'))
TS_mpc = get_many_steps_ts(y_exp, ref_s, step_idx, TOL, verbose, 1, tol_mode);
TS_dat_tmp.ts_s = TS_mpc;
TS_dat_tmp.name = 'MPC with $\mathcal{H}^{-1}$';
TS_dat_cell{end+1} = TS_dat_tmp;

h5 = plot(y_exp.Time, y_exp.Data, '--k');
h5.DisplayName = 'MPC (inv H \& drft)';
% leg = legend([h1, h2, h3,h4, h5]);
leg = legend([h2, h3,h4, h5]);


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

% -- Build up the body of the table. Outer loop is for each row.
% Inner loop is for each experiment (columns).
for k = 2:length(ref_s)
  
  delta_ref = ref_s(k) - ref_s(k-1);
  str_ref_cols = sprintf('&%.2f & %.2f', ref_s(k), delta_ref);
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
    fid = fopen(fullfile(PATHS.MPCJ_root, 'latex', 'manystepsdata.tex'), 'w+');
    fprintf(fid, '%s', S);
    fclose(fid);
end