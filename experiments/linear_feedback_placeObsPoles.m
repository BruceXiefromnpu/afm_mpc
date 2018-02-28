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

addpath('functions')
% PATH_sim_model       = pwd;  % for simulink simulations

% ---- Paths for shuffling data to labview and back. ------
%labview reads data here
controlDataPath = fullfile(PATH_step_exp, controlParamName); 
% labview saves experimental results/data here
dataOut_path    = fullfile(PATH_step_exp, outputDataName); 
% labview reads desired trajectory here
refTrajPath     = fullfile(PATH_step_exp, refTrajName); 
% location of the vi which runs the experiment.
vipath ='C:\Users\arnold\Documents\MATLAB\afm_mpc_journal\labview\play_AFMss_integral_trajTrack_nod.vi';


% ---------- Load Parametric Models  -----------
modFitPath    = 'x-axis_sines_info_out_12-10-2017-04.mat';

load(fullfile(PATH_sysid, modFitPath), 'modelFit')
FitNum = 'sys12_2';
FitNum = 'sys_x_eigen';

TOL = .01;
umax = 5;

SYS = modelFit.(FitNum);
Ts  = SYS.Ts;
%%
SYS.InputDelay = 10;
Nd  = SYS.InputDelay;

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
Rw = modelFit.Rw*5;
Qw = modelFit.Qw;
L1 = dlqr(SYS.A', SYS.C', modelFit.Qw, Rw)';
syscl = ss(SYS.A - L1*SYS.C, SYS.B, SYS.C, 0, Ts);
pzplot(syscl)
C1 = ctrb(SYS);

% sys2 = decimateTF(SYS);
% sys3 = pairDecimated(sys2);
% SYS = ss(sys3);
% C2 = ctrb(SYS);
% 
% Gamma_w = C2*(C1\modelFit.Gamma_w);
% 
% Qw = Gamma_w*modelFit.Rw*Gamma_w';
% L2 = dlqr(SYS.A', SYS.C', Qw, Rw)';
% syscl2 = ss(SYS.A - L2*SYS.C, SYS.B, SYS.C, 0, Ts);
% hold on
% pzplot(syscl2)

%
% 1). Create system with delay but without integral action
PLANT_NOM = absorbDelay(SYS);
PLANT = PLANT_NOM;

%  2). system with integral augmentation.
sys_designK_aug = addIntegral(PLANT);

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
trun = Ts*N1;

ref_f_1 = 1.0; % 1.5 to hit slew rate, 1.4 doesn't  
ref_0 = 0;
t = [0:1:N1]'*Ts;

yref = [0*t + ref_f_1];

ref_traj.signals.values = yref;
ref_traj.time           = t;
rw = 8.508757290909093e-08;
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
gam_s = [1, 1, 1, 1, 1]; % factors to increase cmplx mode freqs by
alp_s = [.85 .5 .4 .4, .4]; % desired damping of cmplx modes
rho_s = [1.2, 1]; % factors to shift real modes by

p_sort = sort_by_w(pole(SYS));
z_sort = sort_by_w(zero(SYS));
p1 = p_sort(1); z1 = z_sort(1);
wp1 = abs(log(p1))/Ts;
wz1 = abs(log(z1))/Ts;

rho_s(1) = .98*wz1/wp1

% fdbk is a class which computes the gain and also stores the data used to
% generate K with. 
Kfdbk = fdbk(sys_designK_aug, 'gams', gam_s, 'pint', p_int, 'alps',...
         alp_s, 'rhos', rho_s, 'doDelay', 1, 'rad', 0.3 );
K_aug = Kfdbk.K;      
Ki = K_aug(1);
K  = K_aug(2:end);
if 1
    sys_cl = ss(sys_designK_aug.a - sys_designK_aug.b*K_aug,...
            sys_designK_aug.b, sys_designK_aug.c, 0, Ts);
    figure(10); clf
    pzplot(PLANT);
    hold on
    pzplot(sys_cl)        
end

% 2). Design FeedForward gains.
%  Use Nbar to cancel integrator pole.
Zi = p_int;
Nbar = -Ki/(Zi-1);
% 3). Design estimator gains, L.
%  No real theoretical justification for this. But it seems to work. 
% L = L2
% L = dlqr(sys_obs.a', sys_obs.c', 55.5*(sys_obs.c'*sys_obs.c), 1)';
% L = dlqr(sys_obs.a', sys_obs.c', blkdiag(Qw, zeros(Nd,Nd)+.000005), Rw)';
L = dlqr(sys_obs.a', sys_obs.c', blkdiag(Qw, zeros(Nd,Nd)+.0000071), Rw)';
% gam_s = [1, 1, 1, 1, 1]*1.5; % factors to increase cmplx mode freqs by
% alp_s = [.57 .7 .7 .7, .2]*1; % desired damping of cmplx modes
% rho_s = [1.2, 1]; % factors to shift real modes by
% Lfdbk = fdbk(sys_obs, 'gams', gam_s, 'alps',...
%          alp_s, 'rhos', rho_s, 'doDelay', 1, 'rad', 0.15 , 'obs', 1);
% L = Lfdbk.K;
Nd = 0;

[uss_0, uss_f, x0, xf, xss] = yss2uss(SYS, ref_f_1, ref_0); 

[~, ~, x0_full] = yss2uss(PLANT, ref_f_1, ref_0);
x0_obs   = x0_full; 

% If we start with a non-zero initial condition, need to preload the
% integrator initial condition. If x0=0, then uss_0=0 and x0_int=0.
x0_int   = uss_0 -ref_0*Nbar+K*x0_full;


sim(fullfile(PATH_sim_models, 'AFM_SS_linearFDBK_obs.slx'));
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

F1 = figure(50);clf
H1 = plot(sim_exp, F1);
% figure(51); plot(thenoise.Data)
%%
%----------------------------------------------------
% Save the controller to .csv file for implementation
clear vi; clear e;
% delay before we start tracking, to let any transients out. Somewhere of a
% debugging setting. 
SettleTicks = 20;  
Iters  = 500



% creat and pack data. Then save it. 
tt = t;
yy = yref;
uKx  = yy*Nbar;

[y_ref, uKx, y_uKx] = pack_uKx_y(uKx, yy, tt);

AllMatrix = packMatrix(sys_obs, L, K*1);
saveControlData(AllMatrix, Ki*1, Nbar, Ns, Nd, ref_f_1, y_uKx, controlDataPath, refTrajPath);

% -----------------------RUN THE Experiment--------------------------------
[e, vi] = setupVI(vipath, 'SettleTicks', SettleTicks, 'Iters', Iters, 'Stop', 0,...
            'Tab Control', 1, 'umax', umax, 'outputData BOTH', dataOut_path,...
            'reference traj path', refTrajPath, 'control_data_path', controlDataPath);
vi.Run
% -------------------------------------------------------------------------
        
        
% Now, read in data, and save to structure, and plot.
AFMdata = csvread(dataOut_path);
[y_exp, u_exp, xhat_exp] = unpackExpData_nod(AFMdata, Ns, Ts);
yy = xhat_exp.Data*sys_obs.c';

expOpts = stepExpOpts(linOpts, 'pstyle', 'g', 'name',  'AFM Stage');

afm_exp = stepExp(y_exp, u_exp, expOpts);
H2 = plot(afm_exp, F1);
subplot(2,1,1)

xk = single(zeros(Ns,1));
xk_plus1 = xk;
y_hat_sim = single(zeros(length(y_exp.Time),1));
for k=1:length(y_exp.Time)
    xk = xk_plus1;
    y_hat_sim(k) = single(sys_obs.c)*xk;
    uk = u_exp.Data(k);
    ye = y_exp.Data(k)-y_hat_sim(k);
    xk_plus1 = single(sys_obs.a)*xk + single(sys_obs.b)*single(uk) +...
        single(L)*single(ye);
    
end

plot(y_exp.Time, y_hat_sim, '--', 'LineWidth', 2)
plot(y_exp.Time, yy, 'k:', 'LineWidth', 2)

figure(200)
plot(y_exp.Time, yy-y_hat_sim)

%%
systems.sys_obs = sys_obs;
systems.PLANT   = PLANT;
systems.sys_designK_aug = sys_designK_aug;
systems.Kfdbk = Kfdbk;
save(fullfile(PATH_step_exp, '09_12_2016_linfdbk_slewproblem_trk2.mat'), 'afm_exp', 'sim_exp',...
      'systems')