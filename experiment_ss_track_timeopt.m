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


% PATH_sim_model       = pwd;  % for simulink simulations

% ---- Paths for shuffling data to labview and back. ------
%labview reads data here
controlDataPath = fullfile(PATHS.step_exp, controlParamName); 
% labview saves experimental results/data here
dataOut_path    = fullfile(PATHS.step_exp, outputDataName); 
% labview reads desired trajectory here
refTrajPath     = fullfile(PATHS.step_exp, refTrajName); 
% location of the vi which runs the experiment.
% vipath ='C:\Users\arnold\Documents\labview\ACC2017_archive\play_AFMss_integral_trajTrack.vi';
vipath ='C:\Users\arnold\Documents\MATLAB\afm_mpc_journal\labview\play_AFMss_integral_trajTrack_nod.vi';


% ---------- Load Parametric Models  -----------
% modFitPath    = 'x-axis_sines_info_out_12-10-2017-04.mat';
% load(fullfile(PATHS.sysid, modFitPath), 'modelFit')
modFitPath    = 'ccta_modelData.mat';
load(fullfile(PATHS.sysid, modFitPath))

FitNum    = 'sys12';


TOL = .01;
umax = 10;

% SYS = modelFit.(FitNum);
SYS = PLANT_init_x;
Ts  = SYS.Ts;
Nd  = SYS.InputDelay;
%%
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
PLANT_NOM = absorbDelay(SYS);
PLANT = PLANT_NOM;

%  2). system with integral augmentation.
do_int = 1
if do_int
    sys_designK_aug = addIntegral(PLANT);
else
    sys_designK_aug = PLANT;
end

% 3). Reduced order system for simulation.
sys_obs            = SYS;
% sys_obs.InputDelay = 0;
Ns  = length(sys_obs.b);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                  Design reference "trajectory"                          %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% We will track two setpoints. Number of samples for each setpoint is 800.
N1    = 800;
trun = Ts*N1;

ref_f_1 = 1.0; % 1.5 to hit slew rate, 1.4 doesn't  
ref_0 = 0;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
% Design control/estimator gains. This gains are what we actually         %
% use in both the simulation and experimentant system.                    %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 1). Design Linar control gain, K.
pstyle = 'b';
p_int = 0.35; % integrator pole location
gam_s = [1.0, 1.0, 1, 1, 1]; % factors to increase cmplx mode freqs by
alp_s = [.85 .7 .7 .7, .7]; % desired damping of cmplx modes
rho_s = [1.2, 1]; % factors to shift real modes by

p_sort = sort_by_w(pole(SYS));
z_sort = sort_by_w(zero(SYS));
p1 = p_sort(1); z1 = z_sort(1);
wp1 = abs(log(p1))/Ts;
wz1 = abs(log(z1))/Ts;

rho_s(1) = wz1/wp1

% fdbk is a class which computes the gain and also stores the data used to
% generate K with. 
Kfdbk = fdbk(sys_designK_aug, 'gams', gam_s, 'pint', p_int, 'alps',...
         alp_s, 'rhos', rho_s, 'doDelay', 1, 'rad', 0.5 );
K_aug = Kfdbk.K;
if do_int
Ki = K_aug(1);
K  = K_aug(2:end);
else
    Ki = 0;
    K = K_aug
end
if verbose >0
    sys_cl = ss(sys_designK_aug.a - sys_designK_aug.b*K_aug,...
            sys_designK_aug.b, sys_designK_aug.c, 0, Ts);
    figure(10)
    pzplot(PLANT);
    hold on
    pzplot(sys_cl)        
end

% 2). Design FeedForward gains.
%  Use Nbar to cancel integrator pole.
Zi = p_int;
% Nbar = -Ki/(Zi-1);


% 3). Design estimator gains, L.
%  No real theoretical justification for this. But it seems to work. 


% sys_obs.InputDelay = Nd;
sys_obs = absorbDelay(sys_obs)
L = dlqr(sys_obs.a', sys_obs.c', 55.5*(sys_obs.c'*sys_obs.c), 1)';
% L = dlqr(sys_obs.a', sys_obs.c', (sys_obs.b*sys_obs.b')*10, 1)';
%

[uss_0, uss_f, x0, xf, xss] = yss2uss(SYS, ref_f_1, ref_0); 

[~, ~, x0_full] = yss2uss(PLANT, ref_f_1, ref_0);
x0_obs   = x0_full; 

% If we start with a non-zero initial condition, need to preload the
% integrator initial condition. If x0=0, then uss_0=0 and x0_int=0.
x0_int   = uss_0;

NoiseP = 0;
if 0
t = [0:1:N1]'*Ts;

yref = [0*t + ref_f_1];

ref_traj.signals.values = yref;
ref_traj.time           = t;
else
    [wp_real_x, wz_real_x] = w_zp_real(SYS);
    rho_1 = wz_real_x(1)/wp_real_x(1);
    g_eject = zpk(exp(-wp_real_x(1)*Ts), exp(-wz_real_x(1)*Ts), 1, Ts);
    g_eject = g_eject/dcgain(g_eject);
    sys_eject = minreal(g_eject*SYS);
    sys_sim = SSTools.deltaUkSys(sys_eject);
    xss_sim = SSTools.getNxNu(sys_sim);

    toBisect = TimeOptBisect(sys_sim, 0.1);
    toBisect.TOL = 1e-8;
    [xx, u] = toBisect.time_opt_bisect(xss_sim*0, xss_sim*ref_f_1)


    figure(1)
    plot(u.Data);
    uu = [u.Data; zeros(500,1)];
    tt = [0:1:length(uu)-1]'*Ts;
    uu = lsim(g_eject, uu, tt);
    lsim(SYS, cumsum(uu), tt);
    [yref, tt, xref] = lsim(absorbDelay(SYS), cumsum(uu), tt);
    uref = cumsum(uu);
    
    ref_traj = timeseries(yref, tt);
    ukx = uref + xref*K';
    uKx = timeseries(ukx, tt);
end
sim(fullfile(PATHS.sim_models, 'AFM_SS_linearFDBK_obs_IMPLEMENT_outer_integral_track.slx'));
y_linear = yplant; 
u_linear = u_full;

linOpts = stepExpOpts('pstyle', '--r', 'TOL', TOL, 'y_ref', ref_f_1,...
                      'controller', Kfdbk, 'name',  'Simulation');

sim_exp = stepExp(y_linear, u_linear, linOpts);


F1 = figure(50);clf;
H1 = plot(sim_exp, F1);


if 1
    sys_cl = ss(sys_obs.a - L*sys_obs.c,...
            sys_obs.b, sys_obs.c, 0, Ts);
    figure(20); clf
    pzplot(PLANT);
    title('observer')
    hold on
    pzplot(sys_cl)        
end



%%
%----------------------------------------------------
% Save the controller to .csv file for implementation
clc
clear vi; clear e;
% delay before we start tracking, to let any transients out. Somewhere of a
% debugging setting. 
SettleTicks = 20;  
Iters  = 450


% creat and pack data. Then save it. 
tt = ref_traj.Time;
yy = ref_traj.Data;
% uKx  = yy*Nbar;
Nbar = 1;
[~, ~, y_uKx] = pack_uKx_y(uKx.Data, yy, tt);

AllMatrix = packMatrix(sys_obs, L, K);
saveControlData(AllMatrix, Ki, Nbar, size(L,1), Nd, ref_f_1, y_uKx, controlDataPath, refTrajPath);

% -----------------------RUN THE Experiment--------------------------------
[e, vi] = setupVI(vipath, 'SettleTicks', SettleTicks, 'Iters', Iters, 'Stop', 0,...
            'Tab Control', 1, 'umax', 4, 'outputData BOTH', dataOut_path,...
            'reference traj path', refTrajPath, 'control_data_path', controlDataPath);
vi.Run
% -------------------------------------------------------------------------
        
%        
% Now, read in data, and save to structure, and plot.
AFMdata = csvread(dataOut_path);
[y_exp, u_exp, xhat_exp] = unpackExpData_nod(AFMdata, Ts);
yy = xhat_exp.Data*sys_obs.c';

expOpts = stepExpOpts(linOpts, 'pstyle', 'g', 'name',  'AFM Stage');

afm_exp = stepExp(y_exp, u_exp, expOpts);
H2 = plot(afm_exp, F1);
subplot(2,1,1)
plot(y_exp.Time, yy, 'k:')




%%
systems.sys_obs = sys_obs;
systems.PLANT   = PLANT;
systems.sys_designK_aug = sys_designK_aug;
systems.Kfdbk = Kfdbk;
save(fullfile(PATH_step_expt, '09_12_2016_linfdbk_slewproblem_trk2.mat'), 'afm_exp', 'sim_exp',...
      'systems')
