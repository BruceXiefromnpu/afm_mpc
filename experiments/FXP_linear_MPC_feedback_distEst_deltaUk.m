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
%%% plants = CanonPlants.plants_with_drift_inv(with_hyst);

% [plants, frf_data] = CanonPlants.plants_drift_inv_hyst_sat();
[plants, frf_data] = CanonPlants.plants_ns14;
% plants.PLANT = plants2.PLANT; %*dcgain(plants.Gvib)/dcgain(plants2.PLANT);
% plants.gdrift = plants2.gdrift;

Ts  = plants.SYS.Ts;
if md == 2
  plants.gdrift = zpk([], [], 1, Ts);
  plants.gdrift_inv = zpk([], [], 1, Ts);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                  Design reference "trajectory"                          %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Get a ref trajectory to track.
N    = 800;
r1 = .65;
r2 = -6;
trajstyle =1;
if trajstyle == 1
  step_ref = StepRef(r1, N);
  yref = step_ref.yref;
elseif trajstyle == 2
  step_ref = StepRef([r1, r2], N);
  yref = step_ref.yref;
elseif trajstyle == 3
  ref_path = fullfile(save_root, 'many_steps_short.mat');
  ref_dat = load(ref_path);
  
  yref = CanonRefTraj.ref_traj_load(ref_path)
elseif trajstyle == 4
  copyfile('many_steps_rand.mat', fullfile(save_root, 'many_steps_rand.mat'));
  ref_path = fullfile(save_root, 'many_steps_rand.mat');
  yref = CanonRefTraj.ref_traj_load(ref_path);
  yref.Data = yref.Data/4;
elseif trajstyle == 5
  step_root = fullfile(PATHS.exp, 'step-exps');
  
  load(fullfile(step_root, 'many_steps_data_rand_ymax7.mat'), 'step_ref');
  yref = step_ref.yref;
end

rw = 8.508757290909093e-07;
rng(1);
thenoise = timeseries(mvnrnd(0, rw, length(yref.Time))*1, yref.Time);


F1 = figure(60); clf
subplot(3,1,1)
hold on, grid on;
step_ref.plot(F1, '-k', 'LineWidth', 0.5)

F11 = figure(61); %clf
step_ref.plot(F11);
step_ref.plot_settle_boundary(F11, TOL, tol_mode);
% legend([h1(1)])



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
du_max_orig = StageParams.du_max;
if md == 1
  du_max = du_max_orig/norm(plants.gdrift_inv, Inf);
else
  du_max = du_max_orig;
end


can_cntrl = CanonCntrlParams_ns14(plants.SYS);
[Q1, R0, S1] = build_control(plants.sys_recyc, can_cntrl);
gam_lin = 3;
gam_mpc = .1;
R1 = R0 + gam_mpc;

K_lqr = dlqr(plants.sys_recyc.a, plants.sys_recyc.b, Q1, R0+gam_lin, S1);
K_lqr2 = dlqr(plants.sys_recyc.a, plants.sys_recyc.b, Q1, R0+gam_mpc, S1);
F = figure(11); clf;
opts.pcolor = 'b';
opts.pstyle = 'x';
opts.zcolor = 'b';
opts.zstyle = 'o';

opts.doZero = 1;
pzplotCL(plants.sys_recyc, K_lqr, [], F, opts);
opts.pcolor = 'r';
opts.pstyle = '*';
opts.zcolor = 'r';
pzplotCL(plants.sys_recyc, K_lqr2, [], F, opts);
%
sys_cl = SSTools.close_loop(plants.sys_recyc, K_lqr);
N_mpc = 12;

Qp = dare(plants.sys_recyc.a, plants.sys_recyc.b, Q1, R1, S1);
nu = 1;
mpcProb = condensedMPCprob_OA(plants.sys_recyc, N_mpc, Q1, Qp, R1, S1);
% mpcProb.lb = zeros(N_mpc,1)-du_max;
% mpcProb.ub = zeros(N_mpc,1)+du_max;
CON = CondenCon([], [], N_mpc);
CON.add_input_con('box', [-du_max, du_max]);
mpcProb.CON = CON;

Hmpc = mpcProb.H; Mmpc = mpcProb.M;
maxIter = 20;
fprintf('condition of H = %.1f\n', mpcProb.kappa);

fgm_fp = FGMprob_1(plants.sys_recyc, N_mpc, Q1, Qp, R1, S1, du_max, maxIter);
I_H_mpc = fgm_fp.I_HL;
ML_x0  = fgm_fp.ML;
beta    = fgm_fp.beta;

% ------------------------- Observer Gain ---------------------------------
can_obs_params = CanonObsParams_01();
[sys_obsDist, L_dist] = build_obs(plants.SYS, can_obs_params);

if 1
  
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
end

% 2). Design FeedForward gains.
[Nx, Nu] = SSTools.getNxNu(plants.sys_recyc);

sims_fpl = SimAFM(plants.PLANT, K_lqr, Nx, sys_obsDist, L_dist, du_max, false, 'thenoise', thenoise);
sims_fpm = SimAFM(plants.PLANT, mpcProb, Nx, sys_obsDist, L_dist, du_max, false, 'thenoise', thenoise);

% if 1
%   sims_fpl.r = plants.hyst.r;
%   sims_fpl.w = plants.hyst.w;
%   sims_fpl.rp = plants.hyst.rp;
%   sims_fpl.wp = plants.hyst.wp;
%   sims_fpl.gdrift_inv = plants.gdrift_inv;
%   sims_fpl.gdrift = plants.gdrift;
% end
if 1
  sims_fpl.r = plants.hyst_sat.r;
  sims_fpl.w = plants.hyst_sat.w;
  sims_fpl.rp = plants.hyst_sat.rp;
  sims_fpl.wp = plants.hyst_sat.wp;
  sims_fpl.gdrift_inv = plants.gdrift_inv;
  sims_fpl.gdrift = plants.gdrift;
  sims_fpl.d = plants.hyst_sat.d;
  sims_fpl.ws = plants.hyst_sat.ws;
  sims_fpl.dp = plants.hyst_sat.dp;
  sims_fpl.wsp = plants.hyst_sat.wsp;
end

[y_lin_fp_sim, U_full_fp_sim, U_nom_fp_sim, dU_fp_sim, Xhat_fp] = sims_fpl.sim(yref);
% ts_lfp = settle_time(y_lin_fp_sim.Time(1:800), y_lin_fp_sim.Data(1:800), r1, TOL*r1);
% fprintf('linear fp settle-time = %.3f [ms]\n', ts_lfp*1000);
% fprintf('perc increase over time-optimal: %.3f\n', (ts_lfp/ts_to)*100);

linOpts = stepExpDuOpts('pstyle', '--k', 'TOL', TOL, 'step_ref', step_ref,...
                      'controller', K_lqr, 'name',  'FP lin Sim.');

sim_exp = stepExpDu(y_lin_fp_sim, U_full_fp_sim, dU_fp_sim, linOpts);
Ts_vec_lfp = sim_exp.settle_time(TOL, tol_mode, 1);
fprintf('Total linear fp settle-time = %.3f [ms]\n', sum(Ts_vec_lfp)*1000);

F1 = figure(60);
h1 = sim_exp.plot(F1, 'umode', 'both');
legend([h1(1)])
h12 = sim_exp.ploty(F11)
legend([h12(1)])
%
figure(70); clf
% du_full = diff(U_full_fp_sim.Data);
% du_full(end+1) = du_full(end);
hold on
% plot(U_full_fp_sim.Time, du_full, '--g')
% plot(dU_fp_sim.Time, dU_fp_sim.Data, 'r')
xlm = xlim();
% plot(xlm, [du_max_orig, du_max_orig], ':k')
% plot(xlm, -[du_max_orig, du_max_orig], ':k')
% legend('du (actual)', 'du (nominal)')
% grid on
F71 = figure(71); clf
% plotState(Xhat_fp, F71);

% -------------------- Setup Fixed stuff -----------------------------

A_obs_cl = sys_obsDist.a - L_dist*sys_obsDist.c;
fprintf('A_cl needs n_int = %d\n', ceil(log2(max(max(abs(A_obs_cl))))) + 1);
fprintf('L needs n_int = %d\n', ceil(log2(max(abs(L_dist)))) + 1);
fprintf('Nx needs n_int = %d\n', ceil(log2(max(abs(Nx)))) + 1);
fprintf('K needs n_int = %d\n', ceil(log2(max(abs(K_lqr)))) + 1);
fprintf('B needs n_int = %d\n', ceil(log2(max(abs(sys_obsDist.b))))+2);

nw = 32;
nf = 26;

du_max_fxp = fi(du_max, 1, 32, 26);
K_fxp = fi(K_lqr, 1, nw,32-10);
Nx_fxp = fi(Nx, 1, 32, 30);
L_fxp = fi(L_dist, 1, 32, 30);

sys_obs_fxp.a = fi(sys_obsDist.a -L_dist*sys_obsDist.c, 1, nw, nw-7);
sys_obs_fxp.b = fi(sys_obsDist.b, 1, nw, 29);
sys_obs_fxp.c = fi(sys_obsDist.c, 1, nw, 28);

% --------------------  Fixed Linear stuff -----------------------------

sims_fxpl = SimAFM(plants.PLANT, K_fxp, Nx_fxp, sys_obs_fxp, L_fxp, du_max_fxp,...
  true, 'nw', nw, 'nf', nf, 'thenoise', thenoise);
% if 1
%   sims_fxpl.r = plants.hyst.r;
%   sims_fxpl.w = plants.hyst.w;
%   sims_fxpl.rp = plants.hyst.rp;
%   sims_fxpl.wp = plants.hyst.wp;
%   sims_fxpl.gdrift_inv = plants.gdrift_inv;
%   sims_fxpl.gdrift = plants.gdrift;
% end
sims_fxpl.r = plants.hyst_sat.r;
sims_fxpl.w = plants.hyst_sat.w;
sims_fxpl.rp = fi(plants.hyst_sat.rp, 1, 16, 11);
sims_fxpl.wp = fi(plants.hyst_sat.wp, 1, 16, 11);
sims_fxpl.d = plants.hyst_sat.d;
sims_fxpl.ws = plants.hyst_sat.ws;
sims_fxpl.dp = fi(plants.hyst_sat.dp, 1, 16, 11);
sims_fxpl.wsp = fi(plants.hyst_sat.wsp, 1, 16, 11);
sims_fxpl.gdrift_inv = plants.gdrift_inv;
sims_fxpl.gdrift = plants.gdrift;


[y_fxpl, U_full_fxpl, U_nom_fxpl, dU_fxpl, Xhat_fxpl] = sims_fxpl.sim(yref);
fxpl_Opts = stepExpDuOpts('pstyle', '-r', 'TOL', TOL, 'step_ref', step_ref,...
                      'controller', K_lqr, 'name',  'FXP lin Sim.');
sim_exp_fxpl = stepExpDu(y_fxpl, U_full_fxpl, dU_fxpl, fxpl_Opts);
Ts_vec_fxpl = sim_exp_fxpl.settle_time(TOL, tol_mode, 1);
fprintf('Total linear fxp settle-time = %.3f [ms]\n', sum(Ts_vec_fxpl)*1000);

h2 = plot(sim_exp_fxpl, F1, 'umode', 'both');
legend([h2(1)])
figure(F11)
h22 = sim_exp_fxpl.ploty(F11);
legend([h12, h22]);

[~, F71] = plotState(Xhat_fxpl, F71, [], [], '--r');
fprintf('max of Xhat = %.2f\n', max(abs(Xhat_fxpl.Data(:))));
fprintf('max of M*Xhat = %.2f\n', max(max(abs(ML_x0*Xhat_fxpl.Data'))));

% Simulate current 
figure(1000); clf
Ipow = lsim(plants.G_uz2powI, U_full_fxpl.Data, U_full_fxpl.Time);
hlin_Ipow = plot(U_full_fxpl.Time, Ipow, '--k');
hlin_Ipow.DisplayName = 'Linear Current';
hold on, grid on;

% ----------------------------------------------------------------------- %
% --------------------- MPC, fgm fixed-point ---------------------------- %
nw_fgm = 32;
nf_fgm = 28;

fgm_fxp = FGMprob_fxp_1(plants.sys_recyc, N_mpc, Q1, Qp, R1, S1, du_max,...
           maxIter, nw_fgm, nf_fgm);
fgm_fxp.x_nw = 32;
fgm_fxp.x_nf = 27;
fprintf('Mx0 needs n_int = %d\n', ceil(log2(double(max(abs(fgm_fxp.ML(:))))))+1);
fprintf('I_HL needs n_int = %d\n', ceil(log2(double(max(abs(fgm_fxp.I_HL(:))))))+1);

sims_fxpm = SimAFM(plants.PLANT, fgm_fxp, Nx_fxp, sys_obs_fxp, L_fxp, du_max_fxp,...
                   true, 'nw', nw, 'nf', nf, 'thenoise', thenoise);
% if 1
%   sims_fxpm.r = plants.hyst.r;
%   sims_fxpm.w = plants.hyst.w;
%   sims_fxpm.rp = plants.hyst.rp;
%   sims_fxpm.wp = plants.hyst.wp;
%   sims_fxpm.gdrift_inv = plants.gdrift_inv;
%   sims_fxpm.gdrift = plants.gdrift;
% end
sims_fxpm.r = plants.hyst_sat.r;
sims_fxpm.w = plants.hyst_sat.w;
sims_fxpm.rp = fi(plants.hyst_sat.rp, 1, 16, 11);
sims_fxpm.wp = fi(plants.hyst_sat.wp, 1, 16, 11);
sims_fxpm.d = plants.hyst_sat.d;
sims_fxpm.ws = plants.hyst_sat.ws;
sims_fxpm.dp = fi(plants.hyst_sat.dp, 1, 16, 11);
sims_fxpm.wsp = fi(plants.hyst_sat.wsp, 1, 16, 11);
sims_fxpm.gdrift_inv = plants.gdrift_inv;
sims_fxpm.gdrift = plants.gdrift;

if 1
  [y_fxpm, U_full_fxpm, U_nom_fxpm, dU_fxpm, Xhat_fxpm, Xerr_fxpm] = sims_fxpm.sim(yref);
  fxpm_Opts = stepExpDuOpts('pstyle', '--g', 'TOL', TOL, 'step_ref', step_ref,...
                          'controller', K_lqr, 'name',  'FXP MPC Simulation');

  sim_exp_fxpm = stepExpDu(y_fxpm, U_full_fxpm, dU_fxpm, fxpm_Opts);
  h3 = plot(sim_exp_fxpm, F1, 'umode', 'both');
  legend([h1(1), h2(1), h3(1)]);
  
  h32 = sim_exp_fxpm.ploty(F11);
  legend([h12, h22, h32]);


  Ts_vec_fxpm = sim_exp_fxpm.settle_time(TOL, tol_mode, 1);
  fprintf('Total MPC FXP settle-time = %.3f [ms]\n', sum(Ts_vec_fxpm)*1000);

figure(1000);
hold on
Ipow = lsim(plants.G_uz2powI, U_full_fxpm.Data, U_full_fxpm.Time);
hmpc_Ipow = plot(U_full_fxpl.Time, Ipow, '--g');
hmpc_Ipow.DisplayName = 'MPC Current';

end


[ts_mat, names] = pretty_print_ts_data(TOL, tol_mode, sim_exp, sim_exp_fxpl, sim_exp_fxpm);

if saveon
  save(fullfile(save_root, 'many_steps_linfxp_sim.mat'), 'sim_exp_fxpl');
       
  save(fullfile(save_root, 'many_steps_mpcfxp_sim.mat'), 'sim_exp_fxpm');
       
end

sims_fxpm.sys_obs_fp = sys_obsDist;
sims_fxpm.sys_obs_fp.a = sys_obsDist.a - L_dist*sys_obsDist.c;

mpc_dat_path = 'Z:\mpc-journal\step-exps\MPCControls01.csv';
traj_path = 'Z:\mpc-journal\step-exps\traj_data.csvtraj_data.csv';
sims_fxpm.write_control_data(mpc_dat_path, yref, traj_path)


return


%%
%--------------------------------------------------------------------------
% -------------- MPC Experiment -------------------------------------------

% Build the u-reset.
if 1
  dry_run = false;
  reset_piezo('t1', 15, 't_final', 25, 'umax', 9, 'k1', 0.55,...
            'verbose', true, 'dry_run', dry_run)
end

%%
SettleTicks = 20000;
% Iters = 2500;
Iters = length(yref.Data)-1;
% create and pack data. Then save it.

[num, den] = tfdata(plants.gdrift_inv);
num = num{1};
den = den{1};

umax = 10.5;
ymax = max(yref.Data)*1.3
clear e;
clear vi;
% -----------------------RUN THE Experiment--------------------------------
vipath =['C:\Users\arnold\Documents\MATLAB\afm_mpc_journal',...
  '\labview\fixed-point-host\play_FXP_AFMss_MPC_distEst_singleAxis.vi'];
if 1
[e, vi] = setupVI(vipath, 'SettleTicks', SettleTicks, 'Iters', Iters,...
   'num', num, 'den', den, 'TF Order', 1*(length(den)-1),...
   'r_s', plants.hyst_sat.rp, 'w_s', plants.hyst_sat.wp, 'N_hyst', 1*length(plants.hyst_sat.rp),...
   'sat_ds', plants.hyst_sat.dp, 'sat_ws', plants.hyst_sat.wsp, 'N_sat', 1*length(plants.hyst_sat.dp),...
   'du_max', du_max,'dry_run', false,...
   'umax', umax, 'ymax', ymax, 'outputDataPath', dataOut_path,...
   'traj_path', traj_path, 'control_data_path', mpc_dat_path);
else
[e, vi] = setupVI(vipath, 'SettleTicks', SettleTicks, 'Iters', Iters,...
          'num', num, 'den', den, 'TF Order', (length(den)-1),...
          'r_s', plants.hyst.rp, 'w_s', plants.hyst.wp,'N_hyst', length(plants.hyst_sat.rp),...
            'umax', umax, 'ymax', ymax, 'outputDataPath', dataOut_path,...
            'dry_run', 0,...
            'traj_path', traj_path, 'control_data_path', mpc_dat_path);
end          
vi.Run

ticks_bench = vi.GetControlValue('Loop Ticks (benchmark)');
fprintf('Actual loop ticks: %d\n', ticks_bench);
% Now, read in data, and save to structure, and plot.
AFMdata = csvread(dataOut_path);

t_exp = (0:size(AFMdata,1)-1)'*Ts;
y_exp = timeseries(AFMdata(:,1), t_exp);
u_exp = timeseries(AFMdata(:, 2), t_exp);
du_exp = timeseries(AFMdata(:,3), t_exp);
ufull_exp = timeseries(AFMdata(:,4), t_exp);

Ipow_exp = timeseries(AFMdata(:,5), t_exp);
xhat_exp = timeseries(AFMdata(:,6:end), t_exp);
yy = xhat_exp.Data*sys_obsDist.c';
expOpts = stepExpDuOpts(linOpts, 'pstyle', '--m', 'name',...
          sprintf('AFM Stage (MPC), gamma=%.1f', gam_mpc);

afm_exp_mpc = stepExpDu(y_exp, ufull_exp, du_exp, expOpts);

Ts_vec_afm_mpc = afm_exp_mpc.settle_time(TOL, tol_mode, 1);
fprintf('Total AFM mpc FXP settle-time = %.3f [ms]\n', sum(Ts_vec_afm_mpc)*1000);


H_mpcexp = plot(afm_exp_mpc, F1, 'umode', 'both');

legend([h1(1), h2(1), h3(1),  H_mpcexp(1) ]);
subplot(3,1,1)
plot(y_exp.Time, yy, ':k')

H_mpcexp2 = afm_exp_mpc.ploty(F11);
legend([h12, h22, h32, H_mpcexp2])

figure(1000); %clf
plot(Ipow_exp.Time, (Ipow_exp.Data/15)*1000)
ylabel('current [mA]')
grid on
title('Current')

[~, F61] = plotState(xhat_exp, F61);
fprintf('Max of experimental Xhat = %.2f\n', max(abs(xhat_exp.data(:))));
[ts_mat, names] = pretty_print_ts_data(TOL, tol_mode, sim_exp, sim_exp_fxpl, sim_exp_fxpm, afm_exp_mpc);
%%

% save('many_steps_data/many_steps_rand_fxpmpc_invHystDrift.mat', 'y_exp', 'u_exp',...
%   'du_exp', 'ufull_exp', 'Ipow_exp', 'yref', 'y_fxpm')
save(fullfile(save_root, 'many_steps_mpc_invHyst_invDrift.mat'), 'afm_exp_mpc');

%%
%--------------------------------------------------------------------------
% --------------------------- LINEAR Experiment ---------------------------

% Build the u-reset.
if 1
  dry_run = false;
  reset_piezo('t1', 15, 't_final', 25, 'umax', 9, 'k1', 0.55,...
            'verbose', true, 'dry_run', dry_run)
end
%
sims_fxpl.sys_obs_fp = sys_obsDist;
sims_fxpl.sys_obs_fp.a = sys_obsDist.a - L_dist*sys_obsDist.c;

fxplin_dat_path = 'Z:\mpc-journal\step-exps\LinControls01.csv';
traj_path = 'Z:\mpc-journal\step-exps\traj_data.csv';
sims_fxpl.write_control_data(fxplin_dat_path, yref, traj_path)

SettleTicks = 20000;
% Iters = length(yref.Data)-1;
Iters = 850;
% create and pack data. Then save it.

[num, den] = tfdata(plants.gdrift_inv);
num = num{1};
den = den{1};

umax = 10.5;
ymax = max(yref.Data)*1.3
clear e;
clear vi;
% -----------------------RUN THE Experiment--------------------------------
vipath =['C:\Users\arnold\Documents\matlab\afm_mpc_journal\',...
  'labview\fixed-point-host\play_FXP_AFMss_LinearDistEst_singleAxis.vi'];
if 1
[e, vi] = setupVI(vipath, 'SettleTicks', SettleTicks, 'Iters', Iters,...
   'num', num, 'den', den, 'TF Order', 1*(length(den)-1),...
   'r_s', plants.hyst_sat.rp, 'w_s', plants.hyst_sat.wp, 'N_hyst', 0*length(plants.hyst_sat.rp),...
   'sat_ds', plants.hyst_sat.dp, 'sat_ws', plants.hyst_sat.wsp, 'N_sat', 0*length(plants.hyst_sat.dp),...
   'du_max', du_max,'dry_run', false,...
   'read_file', true, 'umax', umax, 'ymax', ymax, 'outputDataPath', dataOut_path,...
   'traj_path', traj_path, 'control_data_path', fxplin_dat_path);
else
[e, vi] = setupVI(vipath, 'SettleTicks', SettleTicks, 'Iters', Iters,...
   'num', num, 'den', den, 'TF Order', (length(den)-1),...
   'r_s', plants.hyst.rp, 'w_s', plants.hyst.wp, 'N_hyst', length(plants.hyst.rp),...
   'du_max', du_max,'dry_run', false,...
   'read_file', true, 'umax', umax, 'ymax', ymax, 'outputDataPath', dataOut_path,...
   'traj_path', traj_path, 'control_data_path', fxplin_dat_path);
end


vi.Run

% Now, read in data, and save to structure, and plot.
AFMdata = csvread(dataOut_path);

t_exp = (0:size(AFMdata,1)-1)'*Ts;
y_exp = timeseries(AFMdata(:,1), t_exp);
u_exp = timeseries(AFMdata(:, 2), t_exp);
du_exp = timeseries(AFMdata(:,3), t_exp);
ufull_exp = timeseries(AFMdata(:,4), t_exp);

Ipow_exp = timeseries(AFMdata(:,5), t_exp);
xhat_exp = timeseries(AFMdata(:,6:end), t_exp);
yy = xhat_exp.Data*sys_obsDist.c';
expOpts = stepExpDuOpts(linOpts, 'pstyle', '-b', 'name',  'AFM Stage (Linear)');

afm_exp_lin = stepExpDu(y_exp, ufull_exp, du_exp, expOpts);
Ts_vec_afm_lin = afm_exp_lin.settle_time(TOL, tol_mode, 1);
fprintf('Total AFM lin FXP settle-time = %.3f [ms]\n', sum(Ts_vec_afm_lin)*1000);

H_linexp = plot(afm_exp_lin, F1, 'umode', 'both');
subplot(3,1,1)
% legend([h1(1), h2(1), h3(1), H_mpcexp(1), H_linexp(1)]);
legend([h1(1), h2(1), h3(1), H_linexp(1)]);
plot(y_exp.Time, yy, ':k')

H_linexp2 = afm_exp_lin.ploty(F11);
% legend([h12, h22, h32, H_mpcexp2])
legend([h12, h22, h32, H_linexp2])

figure(1000); clf
plot(Ipow_exp.Time, (Ipow_exp.Data/15)*1000)
ylabel('current [mA]')
grid on
title('Current')

[~, F61] = plotState(xhat_exp, F61);
fprintf('Max of experimental Xhat = %.2f\n', max(abs(xhat_exp.data(:))));
%%

% save('many_steps_data/many_steps_rand_fxpmpc_invHystDrift.mat', 'y_exp', 'u_exp',...
%   'du_exp', 'ufull_exp', 'Ipow_exp', 'yref', 'y_fxpm')
save(fullfile(save_root, 'many_steps_linfxp_invHyst_invDrift.mat'), 'afm_exp_lin')

















