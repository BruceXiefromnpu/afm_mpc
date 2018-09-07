% Implement augmented state space integral control with shift register to
% 'estimate' the delay states.


clear all


% Options
figbase  = 50;
verbose = 0;
controlParamName = 'exp01Controls.csv';
refTrajName      = 'ref_traj_track.csv';
outputDataName = 'exp01outputBOTH.csv';
% Build data paths

addpath('../functions')
addpath('../functions/canon')
addpath('../models')
% PATH_sim_model       = pwd;  % for simulink simulations

% ---- Paths for shuffling data to labview and back. ------
%labview reads data here
controlDataPath = fullfile(PATHS.step_exp, controlParamName);
% labview saves experimental results/data here
dataOut_path    = fullfile(PATHS.step_exp, outputDataName);
% labview reads desired trajectory here
refTrajPath     = fullfile(PATHS.step_exp, refTrajName);
% location of the vi which runs the experiment.

% folder into which to save results. Point
% process_settle_time_data.m here...
experiment_directory = ['many_steps_data_rand_', date, '_01'];
step_exp_root = fullfile(PATHS.exp, 'step-exps');
[status, message ] = mkdir(step_exp_root, experiment_directory);
save_root = fullfile(step_exp_root, experiment_directory);
%%
TOL = .01;
tol_mode = 'abs';
saveon = false;

md = 1;
% ------- Load Plants -----
with_hyst = true;

% [plants, frf_data] = CanonPlants.plants_drift_inv_hyst_sat();
[plants, frf_data] = CanonPlants.plants_ns14;
Ts  = plants.SYS.Ts;
nd = 0
plants.SYS = plants.sys_nodelay;
plants.SYS.InputDelay = nd;
plants.PLANT = plants.SYS;

plants.SYS = absorbDelay(plants.SYS)
plants.sys_recyc = SSTools.deltaUkSys(plants.SYS);
plants.PLANT = absorbDelay(plants.PLANT)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                  Design reference "trajectory"                          %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Get a ref trajectory to track.
N    = 800;
r1 = 1;

step_ref = StepRef(r1, N);
yref = step_ref.yref;
dist_traj = yref;
dist_traj.Data = dist_traj.Data*0 + 1;
yref.Data = yref.Data*0;

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

% Adjust the du_max to account for the gain of gdrift_inv.
du_max = 1000;

can_cntrl = CanonCntrlParams_ns14(plants.SYS);
can_cntrl.rad = .05
can_cntrl.rho_s = [2.5, 1.]
can_cntrl.pint = 0.2
[Q1, R0, S1] = build_control(plants.sys_recyc, can_cntrl);
gam_lin = 100;

K_lqr = dlqr(plants.sys_recyc.a, plants.sys_recyc.b, Q1, R0+gam_lin, S1);

% col_idx = 1;
ColOrd = get(gca,'ColorOrder');
clr = ColOrd(col_idx, :);
col_idx = col_idx+1;
R1 = R0 + gam_mpc;

sys_cl = SSTools.close_loop(plants.sys_recyc, K_lqr);
N_mpc = 12;

% ------------------------- Observer Gain ---------------------------------
can_obs_params = CanonObsParams_01();
[sys_obsDist, L_dist] = build_obs(plants.SYS, can_obs_params);


f10 = figure(10); clf
pzplotCL(plants.sys_recyc, K_lqr, [], f10);

F20 = figure(20); clf
p_plant = pole(plants.PLANT);
z_plant = zero(plants.PLANT);
hpobs = plot(real(p_plant), imag(p_plant), 'xk', 'MarkerSize', 8);
hpobs.DisplayName = 'System Pole';

hold on
hzobs = plot(real(z_plant), imag(z_plant), 'ok', 'MarkerSize', 8);
hzobs.DisplayName = 'System Zero';

opts.pcolor = 'r';
[~, hpcl] = pzplotCL(sys_obsDist, [], L_dist, gcf, opts, 'MarkerSize', 8);
hpcl.DisplayName = 'C.L Obs Pole';

xlim([-0.05, 1.05])
ylim([-0.4, 0.4])

leg1 = legend([hpobs, hzobs, hpcl]);
set(leg1, 'location', 'SouthWest', 'box', 'off', 'FontSize', 14);


% 2). Design FeedForward gains.
[Nx, Nu] = SSTools.getNxNu(plants.sys_recyc);

sims_fpl = SimAFM(plants.PLANT, K_lqr, Nx, sys_obsDist, L_dist, du_max, false, 'thenoise', thenoise);

[y_lin_fp_sim, U_full_fp_sim, U_nom_fp_sim, dU_fp_sim, Xhat_fp] = sims_fpl.sim(yref, dist_traj);

linOpts = stepExpDuOpts('pstyle', 'b', 'TOL', TOL, 'step_ref', step_ref,...
                      'controller', K_lqr, 'name',  'FP lin Sim.');

sim_exp = stepExpDu(y_lin_fp_sim, U_full_fp_sim, dU_fp_sim, linOpts);
Ts_vec_lfp = sim_exp.settle_time(TOL, tol_mode, 1);
fprintf('Total linear fp settle-time = %.3f [ms]\n', sum(Ts_vec_lfp)*1000);

F3 = figure(62); %clf
h12 = plot(y_lin_fp_sim.Time, y_lin_fp_sim.Data, 'Color', clr);
legend([h12(1)])
hold on
%





