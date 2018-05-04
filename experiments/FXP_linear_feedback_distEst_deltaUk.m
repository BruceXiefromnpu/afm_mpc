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
load(fullfile(PATHS.sysid, 'hysteresis/steps_hyst_model.mat'));
load(fullfile(PATHS.sysid, 'FRF_data_current_stage2.mat'))
% r = r;
w = theta_hyst;

umax = 5;

TOL = .01;

%%
% SYS = ss(Gvib);
% PLANT = ss(Gvib);
% PLANT = ss(modelFit.models.G_uz2stage);
% SYS = ss(modelFit.models.G_uz2stage);

PLANT = (ss(modelFit.models.G_uz2stage));
SYS = (ss(modelFit.models.G_uz2stage));



SYS = balreal(SYS);
PLANT = balreal(SYS);
Nx = SSTools.getNxNu(PLANT);
T = diag(1./Nx)/10;
SYS = ss2ss(SYS, T);
PLANT = ss2ss(PLANT, T);

Nd = 9;
SYS.iodelay = 0;
SYS.InputDelay = 0;

PLANT.InputDelay = Nd;
PLANT = absorbDelay(PLANT);

SYS.InputDelay = Nd;
SYS = absorbDelay(SYS);
Ts  = SYS.Ts;

% Nx = SSTools.getNxNu(PLANT);
% T = diag(1./Nx);
% SYS = ss2ss(SYS, T);
% PLANT = ss2ss(PLANT, T);
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
ref_f_1 = 1.5; % 1.5 to hit slew rate, 1.4 doesn't
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

zeta_x = [.85, .7, .4, .4 .4];
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
Nbar = K_lqr*Nx + Nu;
gdrift_inv = 1/gdrift;
sim_struct = struct('K_lqr', K_lqr, 'du_max', du_max, 'PLANT', PLANT,...
             'Nx_recyc', Nx, 'sys_obs', sys_obsDist, 'L', L_dist);
             

[y_linear, U_full, U_nom, dU, Xhat_fp] = sim_AFM(sim_struct, ref_traj);

linOpts = stepExpOpts('pstyle', '-r', 'TOL', TOL, 'y_ref', ref_f_1,...
                      'controller', K_lqr, 'name',  'Simulation');

sim_exp = stepExp(y_linear, U_full, linOpts);

F1 = figure(59); clf
H1 = plot(sim_exp, F1);
subplot(2,1,1)
plot(ref_traj.time, ref_traj.Data, '--k', 'LineWidth', .05);
xlm = xlim();


F2 = figure(60); clf
plot(dU)
hold on
plot(U_nom.Time(1:end-1), diff(U_nom.Data), '--r')
% plot(xlm, [ref_f_2, ref_f_2], ':')

F61 = figure(61); clf
plotState(Xhat_fp, F61);
%
clc

A_obs_cl = sys_obsDist.a - L_dist*sys_obsDist.c;
fprintf('A_cl needs n_int = %d\n', ceil(log2(max(max(abs(A_obs_cl))))) + 1)
fprintf('L needs n_int = %d\n', ceil(log2(max(abs(L_dist)))) + 1)
fprintf('Nx needs n_int = %d\n', ceil(log2(max(abs(Nx)))) + 1)
fprintf('K needs n_int = %d\n', ceil(log2(max(abs(K_lqr)))) + 1)
fprintf('B needs n_int = %d\n', ceil(log2(max(abs(sys_obsDist.b))))+2)

nw = 32;
% nf = 27;

sim_struct.du_max = fi(0.198, 1, 32, 27);
sim_struct.K_lqr = K_lqr;
sim_struct.Nx_recyc = fi(Nx, 1, 32, 30);
sim_struct.L = fi(L_dist, 1, 32, 30);

rmfield(sim_struct, 'sys_obs');
sim_struct.sys_obs.a = fi(sys_obsDist.a -L_dist*sys_obsDist.c, 1, nw, nw-7);
sim_struct.sys_obs.b = fi(sys_obsDist.b, 1, nw, 29);
sim_struct.sys_obs.c = fi(sys_obsDist.c, 1, nw, 28);
sim_struct.x0_obs = sys_obsDist.b*0;
Ident_obs = IDENT_obs;

% sim_struct.du_max = fi(0.198, 1, nw, nf);
sim_struct.K_lqr = fi(K_lqr, 1, nw,32-10);
% sim_struct.Nx_recyc = fi(Nx, 1, nw, nf);
% sim_struct.L = fi(L_dist, 1, nw, nf);

% rmfield(sim_struct, 'sys_obs');
% sim_struct.sys_obs.a = fi(A_obs_cl , 1, nw, nf-2);
% sim_struct.sys_obs.b = fi(sys_obsDist.b, 1, nw, nf);
% sim_struct.sys_obs.c = fi(sys_obsDist.c, 1, nw, nf);
% sim_struct.x0_obs = sim_struct.sys_obs.b*0;
% Ident_obs = fi(IDENT_obs, 1,32,1);
% , eNs_12
C_ydist = eNs_12
x0 = PLANT.b*0;
sim('PREP_FXP_AFMss_obshas_uk')
figure(59)
subplot(2,1,1)
plot(Y, '-b')
subplot(2,1,2)
plot(U, '-b')
plotState(Xhat, F61, [], [], '--')



% AllMatrix2 = [sys_obsDist.b(:), sys_obsDist.c(:), L_dist(:), K_lqr(:),  Nx(:), A_obs_cl, Xhat_fp.Data(20,:)'];
% lv_unittest_path = 'C:\Users\arnold\Documents\MATLAB\afm_mpc_journal\labview\UnitTests\fpga_harnesses';
% test_dataPath = fullfile(lv_unittest_path, 'all_matrix_cols.csv');
% dlmwrite(test_dataPath, AllMatrix2, 'delimiter', ',', 'precision', 12);

%%
%----------------------------------------------------
% Build the u-reset.
% addpath('C:\Users\arnold\Documents\MATLAB\afm_mpc_journal\modelFitting\hysteresis')
% t1 = 10;
% t_final = 20;
% umax = 5;
% k1 = 0.45;
% u_reset = PIHyst.gen_reset_u(t1, t_final, Ts, k1, umax);
% figure(100)
% plot(u_reset)
% grid on
if 1
  reset_piezo();
end
%%
% Save the controller to .csv file for implementation
clear vi; clear e;
% delay before we start tracking, to let any transients out. Somewhere of a
% debugging setting.
SettleTicks = 20000;
Iters = length(yref)-1;

Iters = 600;
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
[y_exp, u_exp, I_exp, xhat_exp] = unpackExpData_nod(AFMdata, Ts);
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
plot(I_exp.Time, (I_exp.Data/15)*1000)
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




















