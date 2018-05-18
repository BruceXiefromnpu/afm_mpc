% Implement augmented state space integral control with shift register to
% 'estimate' the delay states.

% clc
clear
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

% where the different experiments are stored.
save_root = fullfile(PATHS.exp, 'experiments', 'many_steps_data')

% ---------- Load Parametric Models  -----------

load(fullfile(PATHS.sysid, 'hysteresis/steps_hyst_model.mat'));
r = r;
w = theta_hyst;
% load('C:\Users\arnold\Documents\MATLAB\afm_mpc_journal\modelFitting\hysteresis\steps_hyst_model_withsat.mat')

load(fullfile(PATHS.sysid, 'FRF_data_current_stage2.mat'))
whos
umax = 5;

TOL = .01;


%%
SYS = ss(Gvib);
PLANT = ss(Gvib);
% PLANT = ss(modelFit.models.G_uz2stage);
% SYS = ss(modelFit.models.G_uz2stage);

Nd = 9;
SYS.iodelay = 0;
SYS.InputDelay = 0;

PLANT.InputDelay = Nd;
PLANT = absorbDelay(PLANT);

SYS.InputDelay = Nd;
SYS = absorbDelay(SYS);

Ts  = SYS.Ts;



Ns  = length(SYS.b);

%--------------------------------------------------------------------------
%                  Design reference "trajectory"                          %

N    = 800;
r1 = 2;
r2 = -.5;
trajstyle =3;
if trajstyle == 1
  yref = CanonRefTraj.ref_traj_1(r1, N);
elseif trajstyle == 2
    yref = CanonRefTraj.ref_traj_2(r1, r2, N);
elseif trajstyle == 3
  yref = CanonRefTraj.ref_traj_load('many_steps.mat');
elseif trajstyle == 4
  fname_traj = ['C:\Users\arnold\Documents\MATLAB\afm_mpc_journal',...
                '\modelFitting\hysteresis\hyst_input_data_5-4-2018.mat'];
  traj_dat = load(fname_traj);
  t_vec = traj_dat.t_vec;
  u_vec = traj_dat.u_vec;
  yref = timeseries(u_vec*dcgain(PLANT), t_vec);
  t = yref.Time;
  yref = yref.Data;
end

rw = 8.508757290909093e-07;
rng(1);
thenoise = timeseries(mvnrnd(0, rw, length(yref.Time))*0, yref.Time);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
% Design control/estimator gains. These gains are what we actually         %
% use in both the simulation and experimentant system.                    %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% -------------------- Constrained LQR Stuff ------------------------------
du_max   = StageParams.du_max;


can_cntrl = CanonCntrlParams_01(SYS);
[sys_recyc, Q1, R1, S1] = build_control(SYS, can_cntrl);

K_lqr = dlqr(sys_recyc.a, sys_recyc.b, Q1, R1, S1);
sys_cl = SSTools.close_loop(sys_recyc, K_lqr);

% ------------------------- Observer Gain ---------------------------------
can_obs_params = CanonObsParams_01();
[sys_obsDist, L_dist] = build_obs(SYS, can_obs_params);

% ------------------------- Plot Closed loop pz----------------------------
if 1
  f10 = figure(10); clf
  pzplotCL(sys_cl, K_lqr, [], f10);
  figure(20); clf
  pzplot(PLANT);
  title('observer')
  hold on
  opts.pcolor = 'r';
  pzplotCL(sys_obsDist, [], L_dist, gcf, opts);
end

% ------- FeedForward gains.
[Nx, Nu] = SSTools.getNxNu(sys_recyc);
Nbar = K_lqr*Nx + Nu;
gdrift_inv = 1/gdrift;

sims_fp = SimAFM(PLANT, K_lqr, Nx, sys_obsDist, L_dist, du_max, false);
if 1
  sims_fp.r = r;
  sims_fp.w = w;
  sims_fp.rp = rp;
  sims_fp.wp = wp;
  sims_fp.gdrift_inv = gdrift_inv;
  sims_fp.gdrift = gdrift;
end


[y_lin_fp_sim, U_full_fp_sim, U_nom_fp_sim, dU_fp_sim] = sims_fp.sim(yref);

linOpts = stepExpOpts('pstyle', '-r', 'TOL', TOL, 'y_ref', yref.Data(1),...
                      'controller', K_lqr, 'name',  'Simulation');

sim_exp = stepExpDu(y_lin_fp_sim, U_full_fp_sim, dU_fp_sim, linOpts);

if saveon
  save(fullfile(save_root, 'many_steps_linfp_sim.mat'),...
       'y_lin_fp_sim', 'U_full_fp_sim', 'U_nom_fp_sim', 'dU_fp_sim')
end
F1 = figure(56); clf
H1 = plot(sim_exp, F1, 'umode', 'both');
subplot(3,1,1)
plot(yref.time, yref.Data, '--k', 'LineWidth', .05);
xlm = xlim();


F2 = figure(60); clf
plot(dU_fp_sim.Time, dU_fp_sim.Data)
hold on
plot(U_full_fp_sim.Time(1:end-1), diff(U_full_fp_sim.Data), '--r')



return

%%
%----------------------------------------------------
% Build the u-reset.
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
Iters = length(yref)-1;

Iters = min(Iters, length(yref)-1);
Iters = length(yref)-1;

% creat and pack data. Then save it.
tt = t;
yy = yref.Data;
uKx  = yy*Nbar;

[y_ref, uKx, y_uKx] = pack_uKx_y(uKx, yy, tt);
[num, den] = tfdata(gdrift_inv);
num = num{1};
den = den{1};

AllMatrix = packMatrixDistEst(sys_obsDist, L_dist, K_lqr, Nx);
saveControlData(AllMatrix, 0, 0, Ns, Nd, ref_f_1, y_uKx, controlDataPath, refTrajPath);
%
clear e;
clear vi;
% -----------------------RUN THE Experiment--------------------------------
vipath ='C:\Users\arnold\Documents\MATLAB\afm_mpc_journal\labview\play_AFMss_fp_distEst_deltaUk.vi';
md = 2;
if md == 1
[e, vi] = setupVI(vipath, 'SettleTicks', SettleTicks, 'Iters', Iters,...
   'num', num, 'den', den, 'TF Order', (length(den)-1),...
    'r_s', rp, 'w-s', wp,'N_hyst', 7,...
    'N_sat', length(dp), 'd_s', dp, 'wsat_s', wsp, 'du_max', du_max,...
            'umax', 7, 'ymax', 5, 'outputData BOTH', dataOut_path,...
            'reference traj path', refTrajPath, 'control_data_path', controlDataPath);
elseif md == 2
[e, vi] = setupVI(vipath, 'SettleTicks', SettleTicks, 'Iters', Iters,...
   'num', num, 'den', den, 'TF Order', (length(den)-1),...
    'r_s', rp, 'w-s', wp,'N_hyst', 7,...
    'N_sat', 0, 'du_max', du_max,...
            'umax', 7, 'ymax', 5, 'outputData BOTH', dataOut_path,...
            'reference traj path', refTrajPath, 'control_data_path', controlDataPath);
end          
vi.Run
% -------------------------------------------------------------------------
%
% Now, read in data, and save to structure, and plot.
AFMdata = csvread(dataOut_path);
[y_exp, u_exp, I_exp, xhat_exp] = unpackExpData_nod(AFMdata, Ts);
yy = xhat_exp.Data*sys_obsDist.c';
du = u_exp; du.Data = du.Data*0;

expOpts = stepExpOpts(linOpts, 'pstyle', '--g', 'name',  'AFM Stage');

afm_exp = stepExpDu(y_exp, u_exp, du, expOpts);
H2 = plot(afm_exp, F1, 'umode', 'both');
subplot(3,1,1)
plot(y_exp.Time, yy, 'k:')

figure(1000); clf
plot(I_exp.Time, (I_exp.Data/15)*1000)
ylabel('current [mA]')
grid on
title('Current')
%%
if md == 2
save(fullfile(save_root, 'many_steps_hyst_nosat_R4.mat'), 'y_exp', 'u_exp', 'I_exp', 'wp', 'rp')
end

















