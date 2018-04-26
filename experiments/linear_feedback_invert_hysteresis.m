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

% modFitPath = 'x-axis_sines_info_out_2-8-2018-01.mat'
% load(fullfile(PATHS.sysid, modFitPath), 'modelFit')

% FRF_data_current_stage2.mat'))
load(fullfile(PATHS.sysid, 'hysteresis/steps_hyst_model.mat'));
% r = r;
w = theta_hyst;

TOL = .01;
umax = 5;

% SYS = modelFit.(FitNum);
%%
do_integrator = 1;
    

SYS = ss(Gvib);

Nd = 9;
SYS.iodelay = 0;
SYS.InputDelay = 0;
PLANT = ss(Gvib*gdrift);
PLANT.InputDelay = Nd;
PLANT = absorbDelay(PLANT);

%    Nd = SYS.InputDelay;
sys_nodelay = SYS;
SYS.InputDelay = Nd;
SYS = absorbDelay(SYS);

%    PLANT = absorbDelay(SYS);
Ts  = SYS.Ts;
% ws = modelFit.frf.w_s;
% Gfrf = modelFit.frf.G_uz2stage;
%    Gfrf = squeeze(modelFit.frf.G_frf);
% F1 = figure(100); clf
% frfBode(Gfrf, ws/2/pi, F1, 'r', 'hz');
% frfBode(SYS,  ws/2/pi, F1, '--k', 'hz');
% plotPZ_freqs(SYS,F1);

%--------------------------------------------------------------------------
% Build models 
% We will create three systems here:
% 1) sys_designK: the S.S. system with delay absorbed as input delay
%    Use this to create the augmented system, (2).
% 2) sys_designK_aug: same as (1) but we augment with integral state.
%    Use this to design K.
% 3) sys_obs: This is the observer system for simulation. It only includes
%    Ns states, with NO delay, because we are implementing a reduced order
%    observer.


% 1). Create system with delay but without integral action


 
%  2). system with integral augmentation.
if do_integrator
    sys_designK_aug = addIntegral(SYS);
else
    sys_designK_aug = SYS;
end

% 3). Reduced order system for simulation.
sys_obs = absorbDelay(SYS);
%
% sys_obs.InputDelay = 0;
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
ref_traj.signals.values = yref;
ref_traj.time           = t;
rw = 8.508757290909093e-09;
rng(1);
thenoise = timeseries(mvnrnd(0, rw, length(t)), t);
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
% Design control/estimator gains. This gains are what we actually         %
% use in both the simulation and experimentant system.                    %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 1). Design Linar control gain, K.
pstyle = 'b';
p_int = 0.8; % integrator pole location
gam_s = [1, 1, 1, 1, 1, 1]; % factors to increase cmplx mode freqs by
alp_s = [.8 .6 .4 .4, .4, .4]; % desired damping of cmplx modes
rho_s = [1, 1, 1, 1 ]; % factors to shift real modes by

p_sort = sort_by_w(pole(SYS));
z_sort = sort_by_w(zero(SYS));
p1 = p_sort(1); z1 = z_sort(1);
wp1 = abs(log(p1))/Ts;
wz1 = abs(log(z1))/Ts;

rho_s(1) = 1.0*wz1/wp1

% fdbk is a class which computes the gain and also stores the data used to
% generate K with. 
Kfdbk = fdbk(sys_designK_aug, 'gams', gam_s, 'pint', p_int, 'alps',...
         alp_s, 'rhos', rho_s, 'doDelay', 1, 'rad', 0.3);
K_aug = Kfdbk.K;      

% K_aug = dlqr(sys_designK_aug.a, sys_designK_aug.b, sys_designK_aug.c'*sys_designK_aug.c, 1)
if 1
    sys_cl = ss(sys_designK_aug.a - sys_designK_aug.b*K_aug,...
            sys_designK_aug.b, sys_designK_aug.c, 0, Ts);
    figure(10); clf
    hold on
    pzplot(sys_cl)        
end
Lfdbk = fdbk(sys_designK_aug, 'gams', gam_s, 'alps',...
         alp_s, 'rhos', rho_s, 'doDelay', 1, 'rad', 0.3, 'obs',true);
L = Lfdbk.K;
% 3). Design estimator gains, L.
% Qw = sys_obs.b*sys_obs.b'*100;
Qw1 = sys_nodelay.b*sys_nodelay.b'*100;
Qw = blkdiag(Qw1, eye(Nd)*500);
L = dlqr(sys_obs.a', sys_obs.c', Qw, 1)';


% 2). Design FeedForward gains.
[uss_0, uss_f, x0, xf, xss] = yss2uss(SYS, ref_f_1, ref_0); 

[~, ~, x0_full] = yss2uss(PLANT, ref_f_1, ref_0);
[~,~, x0_obs] = yss2uss(sys_obs, ref_f_1, ref_0);
if do_integrator
    Ki = K_aug(1);
    K  = K_aug(2:end);
    %  Use Nbar to cancel integrator pole.
    Zi = p_int;
    Nbar = -Ki/(Zi-1);
    % If we start with a non-zero initial condition, need to preload the
    % integrator initial condition. If x0=0, then uss_0=0 and x0_int=0.
    x0_int   = uss_0 -ref_0*Nbar+K*x0_obs;    
else
    K = K_aug;
   [Nx, Nu] = SSTools.getNxNu(SYS);
   Nbar = K*Nx + Nu;
   
   Ki = 0;
   x0_int = 0;
end
gdrift_inv = 1/gdrift;
sim(fullfile(PATHS.sim_models, 'AFM_SS_linearFDBK_inv_hyst_drift.slx'));
y_linear = yplant; 
u_linear = u_full;

linOpts = stepExpOpts('pstyle', '--r', 'TOL', TOL, 'y_ref', ref_f_1,...
                      'controller', Kfdbk, 'name',  'Simulation');

sim_exp = stepExp(y_linear, u_linear, linOpts);

if 1
    sys_cl = ss(sys_obs.a - L*sys_obs.c,...
            sys_obs.b, sys_obs.c, 0, Ts);
    figure(20); clf
    pzplot(PLANT);
    title('observer')
    hold on
    pzplot(sys_cl)        
end

F1 = figure(59); clf
% subplot(2,1,1)
% plot(xhat.Time, yhat, 'b')
H1 = plot(sim_exp, F1);

%%
%----------------------------------------------------
% Build the u-reset.
addpath('C:\Users\arnold\Documents\MATLAB\afm_mpc_journal\modelFitting\hysteresis')
t1 = 10;
t_final = 20;
umax = 5;
k1 = 0.45;
u_reset = PIHyst.gen_reset_u(t1, t_final, Ts, k1, umax);
figure(100)
plot(u_reset)
grid on
%
slewfname_in = 'hyst_reset_datain.csv';
slewfpath_in = fullfile(pwd, slewfname_in);
slewfname_out = 'hyst_reset_data_out.csv';
slewfpath_out = fullfile(pwd, slewfname_out);

csvwrite(slewfpath_in, u_reset);
clear vi;
vipath_reset = 'C:\Users\arnold\Documents\MATLAB\afm_mpc_journal\labview\play_nPoint_id_slew_OL_FIFO.vi'

[e, vi] = setupVI(vipath_reset, 'Abort', 0,...
        'umax', 9, 'data_out_path', slewfpath_out,...
        'traj_in_path', slewfpath_in, 'TsTicks', 1600);
vi.Run
%%
dat = csvread(slewfpath_out);

u_exp = dat(:,1);
yx_exp = dat(:,2); % - dat(1,2);
t_exp = (0:length(u_exp)-1)'*Ts;

figure(101);
plot(t_exp, u_exp*.68)
hold on
plot(t_exp, yx_exp)
grid on

%%    
% Save the controller to .csv file for implementation
clear vi; clear e;
% delay before we start tracking, to let any transients out. Somewhere of a
% debugging setting. 
SettleTicks = 20000;  
Iters = 1599;
Iters = min(Iters, length(yref)-1);


% creat and pack data. Then save it. 
tt = t;
yy = yref;
uKx  = yy*Nbar;

[y_ref, uKx, y_uKx] = pack_uKx_y(uKx, yy, tt);
[num, den] = tfdata(gdrift_inv);
num = num{1};
den = den{1};

AllMatrix = packMatrix(sys_obs, L, K);
saveControlData(AllMatrix, Ki, Nbar, Ns, Nd, ref_f_1, y_uKx, controlDataPath, refTrajPath);

clear e;
clear vi;
% -----------------------RUN THE Experiment--------------------------------
vipath ='C:\Users\arnold\Documents\MATLAB\afm_mpc_journal\labview\play_AFMss_integral_trajTrack_nod.vi';
[e, vi] = setupVI(vipath, 'SettleTicks', SettleTicks, 'Iters', Iters,...
    'num', num, 'den', den, 'TF Order', length(den)-1,...
    'r_s', rp, 'w-s', wp, 'N_hyst', 7,...
            'Tab Control', 1, 'umax', umax, 'outputData BOTH', dataOut_path,...
            'reference traj path', refTrajPath, 'control_data_path', controlDataPath);
vi.Run
% -------------------------------------------------------------------------
        
%
% Now, read in data, and save to structure, and plot.
AFMdata = csvread(dataOut_path);
[y_exp, u_exp, I_exp, xhat_exp] = unpackExpData_nod(AFMdata, Ts);
yy = xhat_exp.Data*sys_obs.c';

expOpts = stepExpOpts(linOpts, 'pstyle', 'm', 'name',  'AFM Stage');

afm_exp = stepExp(y_exp, u_exp, expOpts);
H2 = plot(afm_exp, F1);
subplot(2,1,1)
plot(y_exp.Time, yy, 'k:')

figure(1000)
plot(I_exp.Time, I_exp.Data/15)
grid on
title('Current')
%%


