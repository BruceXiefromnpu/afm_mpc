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
load(fullfile(PATHS.sysid, 'hysteresis/steps_hyst_model.mat'));
% load('C:\Users\arnold\Documents\MATLAB\afm_mpc_journal\modelFitting\hysteresis\steps_hyst_model_withsat.mat')
load(fullfile(PATHS.sysid, 'FRF_data_current_stage2.mat'))
% r = r;
w = theta_hyst;

umax = 5;

TOL = .01;

%%
md = 3;
if md == 1
  PLANT = ss(modelFit.models.G_uz2stage);
  SYS = ss(modelFit.models.G_uz2stage);
elseif md == 2
  [Gvib, gdrift] = eject_gdrift(modelFit.models.G_uz2stage);
  PLANT = ss(Gvib);
  SYS = ss(Gvib);
elseif md == 3
%   [Gvib, ] = eject_gdrift(modelFit.models.G_uz2stage);
%   PLANT = ss(modelFit.models.G_uz2stage);
  PLANT = ss(Gvib);
  SYS = ss(Gvib);
else
  error('choose mode')
end

% PLANT = (ss(modelFit.models.G_uz2stage));
% SYS = (ss(modelFit.models.G_uz2stage));

SYS = balreal(SYS);
% PLANT = balreal(SYS);
Nx = SSTools.getNxNu(SYS);
T = diag(1./Nx)/10;
SYS = ss2ss(SYS, T);
% PLANT = ss2ss(PLANT, T);

Nd = 9;
SYS.iodelay = 0;
SYS.InputDelay = 0;

PLANT.InputDelay = Nd;
PLANT = absorbDelay(PLANT);

SYS.InputDelay = Nd;
SYS = absorbDelay(SYS);
Ts  = SYS.Ts;


%--------------------------------------------------------------------------
% Build models
sys_designK = SYS;

% 3). Reduced order system for simulation.
sys_obs = absorbDelay(SYS);
Ns  = length(sys_obs.b);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                  Design reference "trajectory"                          %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% We will track one setpoints. Number of samples for each setpoint is 800.
N1    = 800;
ref_f_1 = 1; % 1.5 to hit slew rate, 1.4 doesn't
trajstyle = 2;
if trajstyle == 1
    N2 = 800;
    trun = Ts*(N1 + N2);
    ref_f_2 = -6; % 1.5 to hit slew rate, 1.4 doesn't
    ref_0 = 0;
    t1 = [0:1:N1]'*Ts;
    t2 = [N1+1:1:N1+N2]'*Ts;

    t = [t1; t2];
    yref = [0*t1 + ref_f_1;
            0*t2 + ref_f_2];
    ref_traj = timeseries(yref, t);
elseif trajstyle == 2
    t1 = [0:1:N1]'*Ts;
    trun = Ts*(N1);
    t = t1;
    yref = [0*t1 + ref_f_1];
    ref_0 = 0;
    ref_traj = timeseries(yref, t);
elseif trajstyle == 3
  load('many_steps.mat')
  ref_traj = ref_traj_params.ref_traj;
  ref_traj.Data = ref_traj.Data;
  t = ref_traj.Time;
  yref = ref_traj.Data;
end
rw = 8.508757290909093e-07;
rng(1);
thenoise = timeseries(mvnrnd(0, rw, length(t))*0, t);

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
[wp_real_x, wz_real_x] = w_zp_real(sys_designK);
rho_1 = wz_real_x(1)/wp_real_x(1);
rhos_x = [rho_1, 1., 1];

% zeta_x = [.85, .7, .4, .4 .4];
zeta_x = [.8, .7, .4, .4 .4];
gams_x = [1., 1., 1., 1., 1];

rad = 0.25;
pint = 0.8;
sys_recyc=SSTools.deltaUkSys(SYS);

P_x  = getCharDes(sys_recyc, gams_x, pint, zeta_x, rhos_x, rad);
[Chat, Dhat] = place_zeros(sys_recyc, P_x);
gam = 1;
Q1 = Chat'*Chat;
S1 = Chat'*Dhat;
R1 = Dhat'*Dhat + gam; % Plus gamma.

K_lqr = dlqr(sys_recyc.a, sys_recyc.b, Q1, R1, S1);
sys_cl = SSTools.close_loop(sys_recyc, K_lqr);
N_mpc = 12;
Qp = dare(sys_recyc.a, sys_recyc.b, Q1, R1, S1);
nu = 1;
mpcProb = condensedMPCprob_OA(sys_recyc, N_mpc, Q1, Qp, R1, S1);

% rcond(

Hmpc = mpcProb.H; Mmpc = mpcProb.M;
maxIter = 20;
fprintf('condition of H = %.1f\n', mpcProb.kappa);

fgm_fp = FGMprob_1(sys_recyc, N_mpc, Q1, Qp, R1, S1, du_max, maxIter);
I_H_mpc = fgm_fp.I_HL;
ML_x0  = fgm_fp.ML;
beta    = fgm_fp.beta;
if 1
    f10 = figure(10); clf
    pzplotCL(sys_cl, K_lqr, [], f10);
end

% -------------------------------------------------------------------------
% ------------------------- Observer Gain ---------------------------------
%
sys_obs = absorbDelay(SYS);
p_int_d = .85;
  % state disturbance does not work unless we put the deltaUk state in the
  % observer too. It probably could be made to, but I havnt worked that out.
Qw = sys_obs.b*sys_obs.b'*150;
Lx = dlqr(sys_obs.a', sys_obs.c', Qw, 1)';
[L_dist, sys_obsDist, IDENT_obs, eNs_12] = DistEst.output_dist_est(sys_obs,...
                                             Lx, p_int_d);
[Nx_r, Nx_d, Nu_r, Nu_d] = DistEst.steady_state_gains(sys_obs, sys_obs.b*0, 1);

if 1
    figure(20); clf
    pzplot(PLANT);
    title('observer')
    hold on
    opts.pcolor = 'r';
    pzplotCL(sys_obsDist, [], L_dist, gcf, opts);
end

% 2). Design FeedForward gains.
[Nx, Nu] = SSTools.getNxNu(sys_recyc);
gdrift_inv = 1/gdrift;
sims_fpl = SimAFM(PLANT, K_lqr, Nx, sys_obsDist, L_dist, du_max, false);
if 1
  sims_fpl.r = r;
  sims_fpl.w = w;
%   sims_fpl.rp = rp*0;
%   sims_fpl.wp = wp*0;
  sims_fpl.gdrift_inv = gdrift_inv;
  sims_fpl.gdrift = gdrift;
end

             


[y_linear, U_full, ~, dU, Xhat_fp] = sims_fpl.sim(ref_traj);
% sim_AFM(sim_struct, ref_traj);
ts_lfp = settle_time(y_linear.Time(1:800), y_linear.Data(1:800), ref_f_1, TOL*ref_f_1);
fprintf('linear fp settle-time = %.3f [ms]\n', ts_lfp*1000);
linOpts = stepExpOpts('pstyle', '-r', 'TOL', TOL, 'y_ref', ref_f_1,...
                      'controller', K_lqr, 'name',  'FP Simulation');

sim_exp = stepExpDu(y_linear, U_full, dU, linOpts);

F1 = figure(59); clf
plot(sim_exp, F1, 'umode', 'both');
% H1 = plot(sim_exp, F1);
subplot(3,1,1)
hold on, grid on;
plot(ref_traj.time, ref_traj.Data, '--k', 'LineWidth', .05);


F61 = figure(61); clf
plotState(Xhat_fp, F61);
%%

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

sims_fxpl = SimAFM(PLANT, K_fxp, Nx_fxp, sys_obs_fxp, L_fxp, du_max_fxp,...
  true, 'nw', nw, 'nf', nf);
if 1
  sims_fxpl.r = r;
  sims_fxpl.w = w;
  sims_fxpl.rp = rp;
  sims_fxpl.wp = wp;
%   sims_fxpl.gdrift_inv = gdrift_inv;
sims_fxpl.gdrift_inv = g2;
  sims_fxpl.gdrift = gdrift;
end


[y_fxpl, ~, U_fxpl, dU_fxpl, Xhat_fxpl] = sims_fxpl.sim(ref_traj);
fxpl_Opts = stepExpOpts('pstyle', '--k', 'TOL', TOL, 'y_ref', ref_f_1,...
                      'controller', K_lqr, 'name',  'FXP Simulation');
sim_exp_fxpl = stepExpDu(y_fxpl, U_fxpl, dU_fxpl, fxpl_Opts);
plot(sim_exp_fxpl, F1, 'umode', 'both');


plotState(Xhat_fxpl, F61, [], [], '--');
fprintf('max of Xhat = %.2f\n', max(abs(Xhat_fxpl.Data(:))));
fprintf('max of M*Xhat = %.2f\n', max(max(abs(ML_x0*Xhat_fxpl.Data'))));

%%
nw_fgm = 32;
nf_fgm = 28;
fgm_fp = FGMprob_1(sys_recyc, N_mpc, Q1, Qp, R1, S1, du_max, maxIter);
fgm_fxp = FGMprob_fxp_1(sys_recyc, N_mpc, Q1, Qp, R1, S1, du_max,...
           maxIter, nw_fgm, nf_fgm);
fgm_fxp.x_nw = 32;
fgm_fxp.x_nf = 27;


sims_fxpm = SimAFM(PLANT, fgm_fxp, Nx_fxp, sys_obs_fxp, L_fxp, du_max_fxp,...
  true, 'nw', nw, 'nf', nf);

[y_fxpm, U_fxpm, ~, dU_fxpm, Xhat_fxpm, Xerr_fxpm] = sims_fxpm.sim(ref_traj);
fxpm_Opts = stepExpOpts('pstyle', '--g', 'TOL', TOL, 'y_ref', ref_f_1,...
                      'controller', K_lqr, 'name',  'FXP MPC Simulation');
                   
sim_exp_fxpm = stepExpDu(y_fxpm, U_fxpm, dU_fxpm, fxpm_Opts);
plot(sim_exp_fxpm, F1, 'umode', 'both');

ts_mfxp = settle_time(y_fxpm.Time(1:800), y_fxpm.Data(1:800), ref_f_1, TOL*ref_f_1);
fprintf('mpc fp settle-time = %.3f [ms]\n', ts_mfxp*1000);
%%
clc
sims_fxpm.sys_obs_fp = sys_obsDist;
sims_fxpm.sys_obs_fp.a = sys_obsDist.a - L_dist*sys_obsDist.c;

mpc_dat_path = 'C:\Users\arnold\Documents\MATLAB\afm_mpc_journal\data\MPCControls01.csv';
traj_path = 'C:\Users\arnold\Documents\MATLAB\afm_mpc_journal\data\traj_data.csv';
sims_fxpm.write_control_data(mpc_dat_path, ref_traj, traj_path)


%%
%----------------------------------------------------
% Build the u-reset.
% addpath('C:\Users\arnold\Documents\MATLAB\afm_mpc_journal\modelFitting\hysteresis')
if 1
  reset_piezo();
end
%%
% Save the controller to .csv file for implementation
% delay before we start tracking, to let any transients out. Somewhere of a
% debugging setting.
SettleTicks = 20000;
Iters = 600;
Iters = min(Iters, length(yref)-1);

% create and pack data. Then save it.

[num, den] = tfdata(g2);
num = num{1};
den = den{1};

umax = 3;

clear e;
clear vi;
% -----------------------RUN THE Experiment--------------------------------
vipath ='C:\Users\arnold\Documents\MATLAB\afm_mpc_journal\labview\fixed-point-host\play_FXP_AFMss_MPC_distEst_singleAxis.vi';
% [e, vi] = setupVI(vipath, 'SettleTicks', SettleTicks, 'Iters', Iters,...
%    'num', num, 'den', den, 'TF Order', (length(den)-1),...
%     'r_s', rp, 'w-s', wp,'N_hyst', 7, 'du_max', du_max,...
%             'umax', 7, 'ymax', 5, 'outputData BOTH', dataOut_path,...
%             'reference traj path', refTrajPath, 'control_data_path', controlDataPath);
[e, vi] = setupVI(vipath, 'SettleTicks', SettleTicks, 'Iters', Iters,...
          'num', num, 'den', den, 'TF Order', 0*(length(den)-1),...
          'r_s', rp, 'w-s', wp,'N_hyst', 7, 'dry_run', 0,...
            'umax', umax, 'ymax', 2, 'outputDataPath', dataOut_path,...
            'traj_path', traj_path, 'control_data_path', mpc_dat_path);
vi.Run

% Now, read in data, and save to structure, and plot.
AFMdata = csvread(dataOut_path);
% [y_exp, u_exp, I_exp, xhat_exp] = unpackExpData_nod(AFMdata, Ts);
  
t_exp = (0:size(AFMdata,1)-1)'*Ts;
y_exp      = timeseries(AFMdata(:,1), t_exp);

u_exp      = timeseries(AFMdata(:, 2), t_exp);
du_exp      = timeseries(AFMdata(:,3), t_exp);

% Pull out observer state data.    
xhat_exp = timeseries(AFMdata(:,4:end), t_exp)

yy = xhat_exp.Data*sys_obsDist.c';
expOpts = stepExpOpts(linOpts, 'pstyle', '--m', 'name',  'AFM Stage (MPC)');

afm_exp = stepExpDu(y_exp, u_exp, du_exp, expOpts);
H2 = plot(afm_exp, F1, 'umode', 'both');
subplot(3,1,1)
plot(y_exp.Time, yy, ':k')

figure(1000); clf
plot(du_exp.Time, (du_exp.Data/15)*1000)
ylabel('current [mA]')
grid on
title('Current')

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Save the controller to .csv file for implementation
clear vi; clear e;
% delay before we start tracking, to let any transients out. Somewhere of a
% debugging setting.
SettleTicks = 20000;
Iters = length(yref)-1;

Iters = 60;
Iters = min(Iters, length(yref)-1);


% create and pack data. Then save it.
tt = t;
yy = yref;
uKx  = yy*Nbar;

[y_ref, uKx, y_uKx] = pack_uKx_y(uKx, yy, tt);
[num, den] = tfdata(gdrift_inv);
num = num{1};
den = den{1};

AllMatrix = packMatrixDistEst(sys_obsDist, L_dist, K_lqr, Nx);
saveControlData(AllMatrix, 0, 0, Ns, Nd, ref_f_1, y_uKx, controlDataPath, refTrajPath);

clear e;
clear vi;
% -----------------------RUN THE Experiment--------------------------------
vipath ='C:\Users\arnold\Documents\MATLAB\afm_mpc_journal\labview\play_AFMss_fp_distEst_deltaUk.vi';
% [e, vi] = setupVI(vipath, 'SettleTicks', SettleTicks, 'Iters', Iters,...
%    'num', num, 'den', den, 'TF Order', (length(den)-1),...
%     'r_s', rp, 'w-s', wp,'N_hyst', 7, 'du_max', du_max,...
%             'umax', 7, 'ymax', 5, 'outputData BOTH', dataOut_path,...
%             'reference traj path', refTrajPath, 'control_data_path', controlDataPath);
[e, vi] = setupVI(vipath, 'SettleTicks', SettleTicks, 'Iters', Iters,...
   'TF Order', (length(den)-1)*0, 'N_hyst', 0, 'du_max', du_max,...
            'umax', 7, 'ymax', 5, 'outputData BOTH', dataOut_path,...
            'reference traj path', refTrajPath, 'control_data_path', controlDataPath);
vi.Run
% -------------------------------------------------------------------------
%
% Now, read in data, and save to structure, and plot.
AFMdata = csvread(dataOut_path);
[y_exp, u_exp, du_exp, xhat_exp] = unpackExpData_nod(AFMdata, Ts);
yy = xhat_exp.Data*sys_obsDist.c';

expOpts = stepExpOpts(linOpts, 'pstyle', '--g', 'name',  'AFM Stage');

afm_exp = stepExp(y_exp, u_exp, expOpts);
H2 = plot(afm_exp, F1);
subplot(2,1,1)
plot(y_exp.Time, yy, ':k')
%%
subplot(2,1,1)
plot(y_exp.Time, yy, 'k:')

figure(1000); clf
plot(du_exp.Time, (du_exp.Data/15)*1000)
ylabel('current [mA]')
grid on
title('Current')
%%



A_obs_cl = sys_obsDist.a - L_dist*sys_obsDist.c;
A_obs_cl_vec = [];
for k=1:size(A_obs_cl,1)
  A_obs_cl_vec = [A_obs_cl_vec; A_obs_cl(k, :)'];
end

assert(sum(A_obs_cl_vec == A_obs_cl_vec) == length(A_obs_cl_vec));
AllMatrix = [sys_obsDist.b(:); L_dist; K_lqr(:); A_obs_cl_vec; Nx(:)];

% AllMatrix = packMatrixDistEst(sys_obsDist, L_dist, K_lqr, Nx);
% saveControlData(AllMatrix, 0, 0, Ns, Nd, ref_f_1, y_uKx, controlDataPath, refTrajPath);
%
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




















