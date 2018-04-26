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
% load(fullfile(PATHS.sysid, 'FRF_data_current_stage2.mat'))
% r = r;
w = theta_hyst;

umax = 5;

TOL = .01;

%%
SYS = ss(Gvib);
% SYS = ss(modelFit.models.G_uz2stage);

Nd = 9;
SYS.iodelay = 0;
SYS.InputDelay = 0;
PLANT = ss(Gvib);
% PLANT = ss(modelFit.models.G_uz2stage);
PLANT.InputDelay = Nd;
PLANT = absorbDelay(PLANT);

%    Nd = SYS.InputDelay;
sys_nodelay = SYS;

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
ref_f_1 = .5; % 1.5 to hit slew rate, 1.4 doesn't  
if 1
    N2 = 800;
    trun = Ts*(N1 + N2);
    ref_f_2 = 1.5; % 1.5 to hit slew rate, 1.4 doesn't
    ref_0 = 0;
    t1 = [0:1:N1]'*Ts;
    t2 = [N1+1:1:N1+N2]'*Ts;

    t = [t1; t2];
    yref = [0*t1 + ref_f_1;
            0*t2 + ref_f_2];
else
    t1 = [0:1:N1]'*Ts;
    trun = Ts*(N1);
    t = t1;
    yref = [0*t1 + ref_f_1];
    ref_0 = 0;
end
ref_traj = timeseries(yref, t);

rw = 8.508757290909093e-07;
rng(1);
thenoise = timeseries(mvnrnd(0, rw, length(t)), t);

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

if 1
    f10 = figure(10); clf
    pzplotCL(sys_cl, K_lqr, [], f10);
end

% -------------------------------------------------------------------------
% ------------------------- Observer Gain ---------------------------------
% Qw = sys_obs.b*sys_obs.b'*1000;
% Qw1 = sys_nodelay.b*sys_nodelay.b'*800;
% Qw = blkdiag(Qw1, eye(Nd+1)*100);

if 0
% *****seems to work******
[sys_obsDist, IDENT_obs, eNs_12] = distEstObserver(sys_obs);
Qw = sys_obs.b*sys_obs.b'*550;
L_dist = dlqr(sys_obs.a', sys_obs.c', Qw, 1)';
L_dist(end+1) = .3;
else
% Try doing state disturbance and placing obs and dist poles independantly.
p_int_d = 0.3;
Qw = sys_obs.b*sys_obs.b'*250;
Lx = dlqr(sys_obs.a', sys_obs.c', Qw, 1)';
[L_dist, sys_obsDist, IDENT_obs, eNs_12] = DistEst.state_dist_est(sys_obs, Lx, p_int_d);
end


% seems to work
% % % [sys_obsDist, IDENT_obs, eNs_12]= distEstObserver(sys_obs);
% % % Qw = sys_obsDist.b*sys_obsDist.b'*150;
% % % L_dist = dlqr(sys_obsDist.a', sys_obsDist.c', Qw, 1)';
% % % L_dist(end) = .4;

if 1
    figure(20); clf
    pzplot(PLANT);
    title('observer')
    hold on
    opts.pcolor = 'r';
%     opts.p
    pzplotCL(sys_obsDist, [], L_dist, gcf, opts); 
%     xlim([0.5664    0.9976])
%     ylim([-0.42, 0.42])
end
%
% 2). Design FeedForward gains.
[Nx, Nu] = SSTools.getNxNu(sys_recyc);
Nbar = K_lqr*Nx + Nu;
gdrift_inv = 1/gdrift;
sim_struct = struct('K_lqr', K_lqr, 'du_max', du_max, 'PLANT', PLANT,...
             'mpc_on', 0, 'Nx', Nx, 'thenoise', thenoise,...
             'sys_obs', sys_obsDist, 'L', L_dist, 'state_mode', 3);

[y_linear, u_linear, dU] = sim_AFM(sim_struct, ref_traj);

linOpts = stepExpOpts('pstyle', '--r', 'TOL', TOL, 'y_ref', ref_f_1,...
                      'controller', K_lqr, 'name',  'Simulation');

sim_exp = stepExp(y_linear, u_linear, linOpts);

F1 = figure(59); clf
H1 = plot(sim_exp, F1);


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
if 0
  reset_piezo();
end


%%    
% Save the controller to .csv file for implementation
clear vi; clear e;
% delay before we start tracking, to let any transients out. Somewhere of a
% debugging setting. 
SettleTicks = 20000;
Iters = 1500;
Iters = min(Iters, length(yref)-1);


% creat and pack data. Then save it. 
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
[e, vi] = setupVI(vipath, 'SettleTicks', SettleTicks, 'Iters', Iters,...
   'num', num, 'den', den, 'TF Order', (length(den)-1),...
    'r_s', rp, 'w-s', wp,'N_hyst', 7, 'du_max', du_max,...
            'umax', umax, 'outputData BOTH', dataOut_path,...
            'reference traj path', refTrajPath, 'control_data_path', controlDataPath);
vi.Run
% -------------------------------------------------------------------------

% Now, read in data, and save to structure, and plot.
AFMdata = csvread(dataOut_path);
[y_exp, u_exp, I_exp, xhat_exp] = unpackExpData_nod(AFMdata, Ts);
yy = xhat_exp.Data*sys_obsDist.c';

expOpts = stepExpOpts(linOpts, 'pstyle', 'g', 'name',  'AFM Stage');

afm_exp = stepExp(y_exp, u_exp, expOpts);
H2 = plot(afm_exp, F1);
subplot(2,1,1)
plot(y_exp.Time, yy, 'k:')

figure(1000); clf
plot(I_exp.Time, (I_exp.Data/15)*1000)
ylabel('current [mA]')
grid on
title('Current')
%%


