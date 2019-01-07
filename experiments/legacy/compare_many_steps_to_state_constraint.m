% I think it will be interesting to compare the box constraint on
% delta u to a state constrained problem. What is the difference in
% settling time for the many steps experiment?



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
save_root = fullfile(PATHS.exp, 'experiments', 'many_steps_data');

% ---------- Load Parametric Models  -----------
% load(fullfile(PATHS.sysid, 'hysteresis/steps_hyst_model.mat'));
% load('C:\Users\arnold\Documents\MATLAB\afm_mpc_journal\modelFitting\hysteresis\steps_hyst_model_withsat.mat')
% load(fullfile(PATHS.sysid, 'FRF_data_current_stage2.mat'))
% r = r;
% w = theta_hyst;



TOL = .01;
saveon = true;
%%
md = 1;
with_hyst = true;
plants = CanonPlants.plants_with_drift_inv(with_hyst);

Ts  = plants.SYS.Ts;
if md == 2
  plants.gdrift = zpk([], [], 1, Ts);
  plants.gdrift_inv = zpk([], [], 1, Ts);
end

Imax = StageParams.Imax;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                  Design reference "trajectory"                          %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Get a ref trajectory to track.
N    = 800;
r1 = 5;
r2 = -6;
trajstyle =1;
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

% We don't have an inverse drift model here, but to be fair, use
% the same delta_u_max.
du_max = StageParams.du_max_vib;


can_cntrl = CanonCntrlParams_01(plants.SYS);
[Q1, R0, S1] = build_control(plants.sys_recyc, can_cntrl);
gam_lin = 3;
gam_mpc = .1;
R1 = R0 + gam_mpc;

K_lqr = dlqr(plants.sys_recyc.a, plants.sys_recyc.b, Q1, R0+gam_lin, S1);
sys_cl = SSTools.close_loop(plants.sys_recyc, K_lqr);
N_mpc = 300;

Qp = dare(plants.sys_recyc.a, plants.sys_recyc.b, Q1, R1, S1);
nu = 1;
mpcProb = condensedMPCprob_OA(plants.sys_recyc, N_mpc, Q1, Qp, R1, S1);
mpcProb2 = condensedMPCprob_OA(plants.sys_recyc, N_mpc, Q1, Qp, R1, S1);
% mpcProb.lb = zeros(N_mpc,1)-du_max;
% mpcProb.ub = zeros(N_mpc,1)+du_max;

G_delu2powI = ss(plants.gdrift_inv*plants.G_delu2powI);
CON = CondenCon(G_delu2powI, G_delu2powI.b*0, N_mpc);
CON2 = CondenCon([], [], N_mpc);

CON.add_state_con('box', [0.1]);


% du_max = .198/norm(plants.gdrift_inv);
du_max = 0.099;
CON2.add_input_con('box', [du_max]);
mpcProb.CON = CON;
mpcProb2.CON = CON2;


Nx = SSTools.getNxNu(plants.sys_recyc);

GI = G_delu2powI;

% Load ref data:
dat = load('many_steps_rand_longts.mat');

figure(1); clf
TS1_s = [];
TS2_s = [];

ref_s = [0, 5, -5, 6, -6, 4, -4];
length(dat.ref_traj_params.ref_s) %
for k=2:length(ref_s)
  % r1 = dat.ref_traj_params.ref_s(k-1)*1.5;
  % r2 = dat.ref_traj_params.ref_s(k)*1.5;
  
  r1 = ref_s(k-1);
  r2 = ref_s(k);

  x0 = Nx*r1;
  xerr_0 = x0-Nx*r2;

  mpcProb.CON.x0 = mpcProb.CON.x0*0;
  
  [dU, Xerr] = mpcProb.solve(xerr_0);
  [dU2, Xerr2] = mpcProb2.solve(xerr_0);
  
  t = (0:length(dU)-1)*Ts;
  
  [Y1] = lsim(plants.sys_recyc, dU, t, x0);
  [Y2] = lsim(plants.sys_recyc, dU2, t, x0);
  
  ts1 = settle_time(t, Y1, r2, TOL*abs(r2-r1));
  ts2 = settle_time(t, Y2, r2, TOL*abs(r2-r1));
  fprintf('ts-ipow: %.6f, ts-der-bound: %.6f\n', ts1, ts2)
  TS1_s(k-1) = ts1;
  TS2_s(k-1) = ts2;
  
  derz = zpk([1], [], 1, Ts);
  [Y] = lsim(plants.sys_recyc, dU, t, x0);
  [Y2] = lsim(plants.sys_recyc, dU2, t, x0);

  Ipow = lsim(GI, dU, t);
  Ipow2 = lsim(GI, dU2, t);
  
  subplot(3,1,1)
  plot(t, Y1, '-b')
  hold on, grid on
  plot(t, Y2, '--r')

  subplot(3,1,2)
  plot(t, dU, '-b')
  hold on, grid on
  plot(t, dU2, '--r')

  subplot(3,1,3)
  plot(t, Ipow, '-b')
  hold on;
  plot(t, Ipow2, '--r')

  xlm = xlim;
  plot(xlm, [Imax, Imax], ':k')
  plot(xlm, -[Imax, Imax], ':k')
  drawnow
end



subplot(3,1,1)
ylabel('y')
leg1 = legend('state constrained', '$\Delta u$ input constrained');
set(leg1, 'FontSize', 14)
 
subplot(3,1,2)
ylabel('$\Delta u(k)$')

subplot(3,1,3)
ylabel('current [A]')
xlabel('time [s]')
grid on


return
%%
hyst = load(fullfile(PATHS.sysid, 'hysteresis/steps_hyst_model.mat'));


clc
load(fullfile(PATHS.step_exp, 'many_steps_data', ...
              'many_steps_rand_fxplin_invHystDrift.mat'))
k = 1600;

addpath('../modelFitting/hysteresis')
u_hyst = PIHyst.hyst_play_op(ufull_exp.Data, hyst.r, hyst.w, hyst.w*0);
[II, t]= lsim(GI, u_hyst, ufull_exp.Time);

figure(6); clf
plot(t, II)
hold on

plot(Ipow_exp.Time, (Ipow_exp.Data - mean(Ipow_exp.Data(1:200)))/15.15, '--')