% Implement augmented state space integral control with shift register to
% 'estimate' the delay states.

clc
clear all
%  close all

% Options
figbase  = 50;
verbose = 0;
controlParamName = 'exp01Controls.csv';
refTrajName      = 'ref_traj_track.csv';
outputDataName = 'exp01outputBOTH.csv';
% Build data paths

addpath('../functions')
% PATH_sim_model       = pwd;  % for simulink simulations

% ---- Paths for shuffling data to labview and back. ------
%labview reads data here
controlDataPath = fullfile(PATHS.step_exp, controlParamName);
% labview saves experimental results/data here
dataOut_path    = fullfile(PATHS.step_exp, outputDataName);
% labview reads desired trajectory here
refTrajPath     = fullfile(PATHS.step_exp, refTrajName);
% location of the vi which runs the experiment.

% ---------- Load Parametric Models  -----------
% load(fullfile(PATHS.sysid, 'hysteresis/steps_hyst_model.mat'));
% load(fullfile(PATHS.sysid, 'FRF_data_current_stage2.mat'))
% r = r;
% w = theta_hyst;

umax = 5;

TOL = .01;

%%
close all
md = 1;
% --------------- Load Plants -------------------
with_hyst = true;
plants = CanonPlants.plants_with_drift_inv(with_hyst);


Ts  = plants.SYS.Ts;

if md == 2
  plants.gdrift = zpk([], [], 1, Ts);
  plants.gdrift_inv = zpk([], [], 1, Ts);
end

% 3). Reduced order system for simulation.
% sys_obs = absorbDelay(plants.SYS);
% Ns  = length(sys_obs.b);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                  Design reference "trajectory"                          %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Get a ref trajectory to track.
N    = 800;
r1 = 1;
r2 = -6;
trajstyle =4;
if trajstyle == 1
  yref = CanonRefTraj.ref_traj_1(r1, N);
elseif trajstyle == 2
    yref = CanonRefTraj.ref_traj_2(r1, r2, N);
elseif trajstyle == 3
  yref = CanonRefTraj.ref_traj_load('many_steps.mat');
elseif trajstyle == 4
  yref = CanonRefTraj.ref_traj_load('many_steps_rand_longts.mat');
end
rw = 8.508757290909093e-07;
rng(1);
thenoise = timeseries(mvnrnd(0, rw, length(yref.Time))*0, yref.Time);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
% Design control/estimator gains. This gains are what we actually         %
% use in both the simulation and experimentant system.                    %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 3). LQR generation gain.
% -------------------------------------------------------------------------
% -------------------- Constrained LQR Stuff ------------------------------
du_max   = StageParams.du_max;

% Pull out open-loop pole-zero information.
can_cntrl = CanonCntrlParams_01(plants.SYS);
[Q1, R0, S1] = build_control(plants.sys_recyc, can_cntrl);
gam_lin = 3;
gam_mpc = 0.5;
R1 = R0 + gam_mpc;

K_lqr = dlqr(plants.sys_recyc.a, plants.sys_recyc.b, Q1, R0+gam_lin, S1);
sys_cl = SSTools.close_loop(plants.sys_recyc, K_lqr);
if 1
    f10 = figure(10); clf
    pzplotCL(sys_cl, K_lqr, [], f10);
end

% -------------------------------------------------------------------------
% ------------------------- Observer Gain ---------------------------------
%
% ------------------------- Observer Gain ---------------------------------
can_obs_params = CanonObsParams_01();
[sys_obsDist, L_dist] = build_obs(plants.SYS, can_obs_params);
if 1
    figure(20); clf
    pzplot(plants.PLANT);
    title('observer')
    hold on
    opts.pcolor = 'r';
    pzplotCL(sys_obsDist, [], L_dist, gcf, opts);
end

% 2). Design FeedForward gains.
[Nx, Nu] = SSTools.getNxNu(plants.sys_recyc);

isfxp = false;
sims_fpl = SimAFM(plants.PLANT, K_lqr, Nx, sys_obsDist, L_dist, du_max, false);
if 1
  sims_fpl.r = plants.hyst.r;
  sims_fpl.w = plants.hyst.w;
  sims_fpl.rp = plants.hyst.rp;
  sims_fpl.wp = plants.hyst.wp;
  sims_fpl.gdrift_inv = plants.gdrift_inv;
  sims_fpl.gdrift = plants.gdrift;
end

[y_lin_fp_sim, U_full_fp_sim, U_nom_fp_sim, dU_fp_sim, Xhat_fp] = sims_fpl.sim(yref);

linOpts = stepExpOpts('pstyle', '-r', 'TOL', TOL, 'y_ref', yref.Data(1),...
                      'controller', K_lqr, 'name',  'Simulation');

sim_exp = stepExpDu(y_lin_fp_sim, U_full_fp_sim, dU_fp_sim, linOpts);

F1 = figure(59); clf
h1 = plot(sim_exp, F1);
subplot(3,1,1)
plot(yref.time, yref.Data, '--k', 'LineWidth', .05);
xlm = xlim();

F61 = figure(61); clf
plotState(Xhat_fp, F61);

% -------------------- Setup Fixed stuff -----------------------------

A_obs_cl = sys_obsDist.a - L_dist*sys_obsDist.c;
fprintf('A_cl needs n_int = %d\n', ceil(log2(max(max(abs(A_obs_cl))))) + 1)
fprintf('L needs n_int = %d\n', ceil(log2(max(abs(L_dist)))) + 1)
fprintf('Nx needs n_int = %d\n', ceil(log2(max(abs(Nx)))) + 1)
fprintf('K needs n_int = %d\n', ceil(log2(max(abs(K_lqr)))) + 1)
fprintf('B needs n_int = %d\n', ceil(log2(max(abs(sys_obsDist.b))))+2)

nw = 32;
nf = 26;

du_max_fxp = fi(0.198, 1, 32, 27);
K_fxp = fi(K_lqr, 1, nw,32-10);
Nx_fxp = fi(Nx, 1, 32, 30);
L_fxp = fi(L_dist, 1, 32, 30);

sys_obs_fxp.a = fi(sys_obsDist.a -L_dist*sys_obsDist.c, 1, nw, nw-7);
sys_obs_fxp.b = fi(sys_obsDist.b, 1, nw, 29);
sys_obs_fxp.c = fi(sys_obsDist.c, 1, nw, 28);

% --------------------  Fixed Linear stuff -----------------------------
%
clc
sims_fxpl = SimAFM(plants.PLANT, K_fxp, Nx_fxp, sys_obs_fxp, L_fxp, du_max_fxp,...
  true, 'nw', nw, 'nf', nf);

sims_fxpl.r = plants.hyst.r;
sims_fxpl.w = plants.hyst.w;
sims_fxpl.rp = fi(plants.hyst.rp, 1, 16, 11);
sims_fxpl.wp = fi(plants.hyst.wp, 1, 16, 11);
sims_fxpl.gdrift_inv = plants.gdrift_inv;
sims_fxpl.gdrift = plants.gdrift;



[y_fxpl, U_full_fxpl, U_nom_fxpl, dU_fxpl, Xhat_fxpl] = sims_fxpl.sim(yref);
fxpl_Opts = stepExpOpts('pstyle', '--k', 'TOL', TOL, 'y_ref', r1,...
                      'controller', K_lqr, 'name',  'FXP lin Sim.');
sim_exp_fxpl = stepExpDu(y_fxpl, U_full_fxpl, dU_fxpl, fxpl_Opts);

h2 = plot(sim_exp_fxpl, F1, 'umode', 'both');
legend([h1(1), h2(1)])
%
fprintf('Max of U_nom = %.2f\n', max(U_nom_fxpl.Data));
fprintf('Max of U_full = %.2f\n', max(U_full_fxpl.Data));

[~, F61] = plotState(Xhat_fxpl, F61, [], [], '--');
fprintf('max of Xhat = %.2f\n', max(abs(Xhat_fxpl.Data(:))));





sims_fxpl.sys_obs_fp = sys_obsDist;
sims_fxpl.sys_obs_fp.a = sys_obsDist.a - L_dist*sys_obsDist.c;

fxplin_dat_path = 'Z:\mpc-journal\step-exps\FXP_lin_Controls01.csv';
traj_path = 'Z:\mpc-journal\step-exps\traj_data.csv';
sims_fxpl.write_control_data(fxplin_dat_path, yref, traj_path)
%----------------------------------------------------
% Build the u-reset.
%%
if 1
  dry_run = false;
  reset_piezo('t1', 15, 't_final', 25, 'umax', 9, 'k1', 0.55,...
            'verbose', true, 'dry_run', dry_run)
end
%%
% Save the controller to .csv file for implementation
clear vi; clear e;
% delay before we start tracking, to let any transients out. Somewhere of a
% debugging setting.
SettleTicks = 20000;
Iters = length(yref.Data)-1;

% Iters =700;
Iters = min(Iters, length(yref.Data)-1);


clc
[num, den] = tfdata(plants.gdrift_inv);
num = num{1};
den = den{1};

umax = 10;
ymax = max(yref.Data)*1.3
clear e;
clear vi;
% -----------------------RUN THE Experiment--------------------------------
vipath =['C:\Users\arnold\Documents\matlab\afm_mpc_journal\',...
  'labview\fixed-point-host\play_FXP_AFMss_LinearDistEst_singleAxis.vi'];


% [e, vi] = setupVI(vipath, 'SettleTicks', SettleTicks, 'Iters', Iters,...
%    'num', num, 'den', den, 'TF Order', (length(den)-1),...
%     'r_s', rp, 'w-s', wp,'N_hyst', 7, 'du_max', du_max,...
%             'umax', 7, 'ymax', 5, 'outputData BOTH', dataOut_path,...
%             'reference traj path', refTrajPath, 'control_data_path', controlDataPath);
% rp = 0; wp = 0;
[e, vi] = setupVI(vipath, 'SettleTicks', SettleTicks, 'Iters', Iters,...
   'num', num, 'den', den, 'TF Order', (length(den)-1),...
   'r_s', plants.hyst.rp, 'w_s', plants.hyst.wp, 'N_hyst', length(plants.hyst.rp),...
   'du_max', du_max,'dry_run', false,...
   'read_file', true, 'umax', umax, 'ymax', ymax, 'outputDataPath', dataOut_path,...
   'traj_path', traj_path, 'control_data_path', fxplin_dat_path);

vi.Run
% -------------------------------------------------------------------------
%
% Now, read in data, and save to structure, and plot.
% AFMdata = csvread(dataOut_path);
AFMdata = vi.GetControlValue('result_data');

t_exp = (0:size(AFMdata,1)-1)'*Ts;
y_exp = timeseries(AFMdata(:,1), t_exp);
u_exp = timeseries(AFMdata(:, 2), t_exp);
du_exp = timeseries(AFMdata(:,3), t_exp);
ufull_exp = timeseries(AFMdata(:,4), t_exp);

Ipow_exp = timeseries(AFMdata(:,5), t_exp);
% xhat_exp = timeseries(AFMdata(:,6:end), t_exp);
% yy = xhat_exp.Data*sys_obsDist.c';

expOpts = stepExpOpts(linOpts, 'pstyle', '--g', 'name',  'AFM Stage');

afm_exp = stepExpDu(y_exp, ufull_exp, du_exp, expOpts);
H2 = plot(afm_exp, F1);
subplot(3,1,1)
% plot(y_exp.Time, yy, ':k')
subplot(3,1,2)

plot(u_exp.Time, u_exp.data, '--m')

figure(1000); clf
plot(Ipow_exp.Time, (Ipow_exp.Data/15)*1000)
ylabel('current [mA]')
grid on
title('Current')
%%

save('many_steps_data/many_steps_rand_fxplin_invHystDrift.mat', 'y_exp', 'u_exp',...
  'du_exp', 'ufull_exp', 'Ipow_exp', 'yref', 'y_lin_fp_sim')



%%



A_obs_cl = sys_obsDist.a - L_dist*sys_obsDist.c;
A_obs_cl_vec = [];
for k=1:size(A_obs_cl,1)
  A_obs_cl_vec = [A_obs_cl_vec; A_obs_cl(k, :)'];
end

assert(sum(A_obs_cl_vec == A_obs_cl_vec) == length(A_obs_cl_vec));
AllMatrix = [sys_obsDist.b(:); L_dist; K_lqr(:); A_obs_cl_vec; Nx(:)];

clear e;
clear vi;
ns = length(K_lqr)-1;


% test_dataPath = fullfile(lv_unittest_path, 'all_matrix_cols.csv');
dlmwrite(controlDataPath, AllMatrix, 'delimiter', ',', 'precision', 12)


Iters = 500
SettleTicks = 100;
% -----------------------RUN THE Experiment--------------------------------
vipath ='C:\Users\arnold\Documents\MATLAB\afm_mpc_journal\labview\fixed-point-host\play_FXP_AFMss_LinearDistEst_singleAxis.vi';
[e, vi] = setupVI(vipath, 'SettleTicks', SettleTicks, 'Iters', Iters,...
           'umax', 3, 'ymax', 5, 'du_max', du_max, 'Ns', ns, 'x_ref', ref_f_1,...
           'All_matrix', AllMatrix, 'dry_run', false, 'read_write_file', false,...
           'dry_run', false, 'All_matrix', AllMatrix);
            
vi.Run

dat = vi.GetControlValue('result_data');
size(dat)


y_exp2 = dat(:,1);
u_exp2 = dat(:,2);
t = (0:length(y_exp2)-1)'*Ts;
x_exp2 = timeseries(dat(:,4:end), t);
yhat = x_exp2.Data*sys_obsDist.c';
figure(59)
subplot(2,1,1)
h1 = plot(t, y_exp2, '--m')
plot(t, yhat, ':b')
subplot(2,1,2)
h1 = plot(t, u_exp2, '--m')

% plotState(x_exp2);




















