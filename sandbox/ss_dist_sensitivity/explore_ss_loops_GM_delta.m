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

% folder into which to save results. Point
% process_settle_time_data.m here...
experiment_directory = ['many_steps_data_rand_', date, '_01'];
step_exp_root = fullfile(PATHS.exp, 'step-exps');
[status, message ] = mkdir(step_exp_root, experiment_directory);
save_root = fullfile(step_exp_root, experiment_directory);
%%
TOL = .01;
tol_mode = 'abs';
% which simulations to run
do_simlinfp = true;

plotstate = false;
plotpoles = false;
saveon = false;


% ------- Load Plants -----

[plants, frf_data] = CanonPlants.plants_drift_inv_hyst_sat(0);

Ts  = plants.SYS.Ts;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                  Design reference "trajectory"                          %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Get a ref trajectory to track.
N  = 800;
r1 = 5;
r2 = -0;

step_ref = StepRef(r1, N);
yref = step_ref.yref;
dist_traj = yref;
dist_traj.Data = dist_traj.Data*0 + 0;
yref.Data = yref.Data*1;

rw = 2e-07;
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
du_max_orig = StageParams.du_max;
G_recyc = plants.sys_recyc;
G = plants.SYS;

% can_cntrl = CanonCntrlParams_ns14(plants.SYS);
% can_cntrl = CanonCntrlParams_01(G);
% % can_cntrl = can_cntrl.aggressive_params();
% [Q1, R0, S1, P_x] = build_control(G, can_cntrl);
% %
% gam_lin = .1;
% R1 = R0 + gam_lin;
p_int = 0.65; cmplx_rad = 0.85; rho_s = [1.5, 1]; rad = 0.5;
Px = getCharDes_const_sig(G_recyc, p_int, cmplx_rad, rho_s, rad).';
% Px = getCharDes_const_sig(G, p_int, cmplx_rad, rho_s, rad).'
z = tzero(G_recyc);
z = sort_by_w(z(imag(z)~=0))
Px(end-1:end) = z(end-1:end);
% Px(1) = []
[Chat, Dhat] = place_zeros(G_recyc, Px)
% [Chat, Dhat] = place_zeros(G, Px)
% Qw = Chat'*Chat;
Q1 = Chat'*Chat;
% Q1 = G_recyc.c'*G_recyc.c;
S1 = Chat'*Dhat;
R0 = Dhat'*Dhat;
gam_lin = 2000;
R1 = 1;
% R1 = R0+gam_lin;

[K_lqr, Pz] = dlqr(G_recyc.a, G_recyc.b, Q1, R1);

Qp = dare(G_recyc.a, G_recyc.b, Q1, R1);

% Estimator
Qw = G.b*G.b';
% Rw = 1*beta;
Lx = G.a*dlqr(G.a', G.c', Qw, 0.001)';
[Sens, Hyd, Hyr, Hyeta, Loop] = ss_loops_delta(plants.SYS, G_recyc, K_lqr, Lx);

figure(1)
step(Hyr)
%

g_direct = ss(G_recyc.a, G_recyc.b, K_lqr, 0, Ts);
figure(1)
rlocus(Loop)

gamm_ = sqrt(R1/(G_recyc.b'*Pz*G_recyc.b+R1));

Gm_mns = 20*log10(1/(1+gamm_));
Gm_pls = 20*log10(1/(1-gamm_));
fprintf('\n------------------------\n');
fprintf('GM bounds (possibly better)\nUpside GM: %f [dB]\nDownside Gm %f [dB]\n', Gm_pls, Gm_mns);

rad = gamm_;


[~, p] = tfdata(G_recyc, 'v');
% p0 = p(end)
% [~, r] = tfdata(g_direct, 'v');
% r0 = r(end)
% rad = r0/p0
t = 0.001:0.01:2*pi;
x = sin(t);
y = cos(t);

figure(2); clf
nyquist(g_direct, 'b');
hold on
plot(rad*x-1, rad*y, 'r')

plot(x-1, y, '--k')
grid
xlim([-2, 0])
ylim([-1.5, 1.5])

nyquist(Loop, '--b')

[gm, pm] = margin(Loop);
fprintf(['Margins with observer:\nGM:%f [dB]\nPM: %f [deg]\n'], 20*log10(gm), pm);

figure(10)
bode(Loop)
hold on

%%
figure(2)
bode(Sens)
grid on
title('Sens')

figure(3)
rlocus(Loop)

figure(4)
bode(Loop)
grid on
title('Loop')

figure(5)
bode(Hyr)
title('Hyr')
grid on
% ------------------------- Observer Gain ---------------------------------
%%
clc
betas = logspace(log10(1), log10(1e-12), 25);
for beta = betas
Qw = G_recyc.b*G_recyc.b';
Rw = 1*beta;
%%
Lx = dlqr(G_recyc.a', G_recyc.c', Qw, 1)';

% Lx = g.a^3*g.b/(g.c*(g.a^2)*g.b);


[Sens, Hyd, Hyr, Hyeta, Loop] = ss_loops(G_recyc, K_lqr, Lx);

figure(24)
bode(Loop)
hold on, grid on

F_clbode = figure(25); hold on, grid on
bode(Sens)



[Gm, Pm] = margin(Loop);
fprintf('Gm=%f, beta=%.f\n', Gm, beta);

drawnow()
end
%%
if plotpoles
  
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


% ------------------------ Linear floating point simulation ---------------
if do_simlinfp
  sims_fpl = SimAFM(plants.PLANT, K_lqr, Nx, sys_obsDist, L_dist, du_max, false, 'thenoise', thenoise);
  
  [y_lin_fp_sim, U_full_fp_sim, U_nom_fp_sim, dU_fp_sim, Xhat_fp] = sims_fpl.sim(yref, dist_traj);
  
  linOpts = stepExpDuOpts('pstyle', '-b', 'TOL', TOL, 'step_ref', step_ref,...
    'controller', K_lqr, 'name',  'FP lin Sim.');
  
  sim_exp_fpl = stepExpDu(y_lin_fp_sim, U_full_fp_sim, dU_fp_sim, linOpts);
  
  
  Ts_vec_lfp = sim_exp_fpl.settle_time(TOL, tol_mode, 1);
  fprintf('Total linear fp settle-time = %.3f [ms]\n', sum(Ts_vec_lfp)*1000);
  
  h1 = sim_exp_fpl.plot(F_yudu, 'umode', 'both');
  legend([h1(1)])
  
  h12 = sim_exp_fpl.ploty(F_y);
  % y_hry = lsim(Hyr, yref.Data, yref.Time);
  y_hry = lsim(Hyeta, thenoise.Data, thenoise.Time);
  % plot(U_nom_fp_sim.Time, y_hry, '--k')
  legend([h12(1)])
  
  Ry_est = cov(y_lin_fp_sim.Data);
  Pcl = dlyap(Hyeta.a, Hyeta.b*rw*Hyeta.b');
  Ry_calc = Hyeta.c*Pcl*Hyeta.c'; %+ rw
  
  sixsig = 6*sqrt(Ry_calc);
  needed_14volt_cov = 14/512;
  fprintf('six-sigma cov: %.4g\n14v x 512 pix: %.4g\n', sixsig, needed_14volt_cov);
end


if do_simmpcfp
  sims_fpm = SimAFM(plants.PLANT, mpcProb, Nx, sys_obsDist, L_dist, du_max, false, 'thenoise', thenoise);

  [y_mpc_fp_sim, U_full_fpm_sim, U_nom_fpm_sim, dU_fpm_sim, Xhat_fpm] = sims_fpm.sim(yref, dist_traj);
  mpcfpOpts = stepExpDuOpts('pstyle', '-r', 'TOL', TOL, 'step_ref', step_ref,...
    'controller', mpcProb, 'name',  'FP mpc Sim.');
  sim_exp_fpm = stepExpDu(y_mpc_fp_sim, U_full_fpm_sim, dU_fpm_sim, mpcfpOpts);
  
  h1mpc = sim_exp_fpm.plot(F_yudu, 'umode', 'both');
  legend([h1(1), h1mpc(1)])
  
  h12mpc = sim_exp_fpm.ploty(F_y);
  legend([h12(1), h12mpc(1)])

end

if plotstate
  F_state = figure(71); clf
  plotState(Xhat_fp, F_state);
end
















