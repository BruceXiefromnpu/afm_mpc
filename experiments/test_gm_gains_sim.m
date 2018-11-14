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

addpath(fullfile(getMatPath(), 'afm_mpc_journal', 'functions'));
addpath(fullfile(getMatPath(), 'afm_mpc_journal', 'functions' , 'canon'));
addpath(fullfile(getMatPath(), 'afm_mpc_journal', '/models'));
% PATH_sim_model       = pwd;  % for simulink simulations

% ---- Paths for shuffling data to labview and back. ------
%labview reads data here
controlDataPath = fullfile(PATHS.step_exp, controlParamName);
% labview saves experimental results/data here
data_out_path    = fullfile(PATHS.step_exp, outputDataName);
% labview reads desired trajectory here
refTrajPath     = fullfile(PATHS.step_exp, refTrajName);
% location of the vi which runs the experiment.



%%
fprintf('\n\n\n\n')
TOL = 14/512; % max volts by pixels
% TOL = .01;
tol_mode = 'abs';
% which simulations to run

do_sim_linfxp = true;
do_sim_mpcfxp = true;

do_sim_hyst = false;
% do_inv_hyst = false;
% do_drift = false;
% do_invdrift = false;


plotpoles = false;
plotdriftpoles = false;
plotdriftbode = false;
plot_powspec = true;
saveon = true;

% We need four modes in this script:
% 1). Constant sigma, minimum-gamma 
% 2). Constant sigma, rob-optimal
% 3). Choose zeta, minimum-gamma
% 4). Choose zeta, rob-optimal
%
% For all cases, drift and hysteresis compensation remains
% unchanged. 
md = 1;

% !!!!!! These gamas should match the data in table II !!!!!!
if md == 1
  % 1). Constant sigma, minimum-gamma 
  gam_lin = 7.5;
  exp_id_str = 'const-sig-min-gam';
elseif md ==2
% 3). Choose zeta, minimum-gamma
  gam_lin = 3.5;
  gam_mpc = 0.00001;
  exp_id_str = 'choose-zet-min-gam';
elseif md == 3
  gam_mpc = 0.001;
  gam_lin = .0001;
  exp_id_str = 'const-sig-same-sig';
end


% ------- Load Plants -----
with_hyst = true;

[plants, frf_data] = CanonPlants.plants_ns14(9, 2);

Ts  = plants.SYS.Ts;

if plotdriftpoles
  figure(91); clf
  pzplot(plants.gdrift)
  hold on
  pzplot(plants.gdrift_inv);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                  Design reference "trajectory"                          %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Get a ref trajectory to track.
exp_num = '06';
N  = 800;
r1 = 1.705;
r2 = -7;
trajstyle =4;
if trajstyle == 1
  step_ref = StepRef([r1], N);
  yref = step_ref.yref;
  yref.Data = yref.Data*1;
  step_descr = 'single_step';
elseif trajstyle == 2
  step_ref = StepRef([r1, r2], N);
  yref = step_ref.yref;
  yref.Data = yref.Data*1;
  step_descr = 'two_step';  
elseif trajstyle == 4
  step_root = fullfile(PATHS.exp, 'step-exps');
  
  load(fullfile(step_root, 'many_steps_data_rand_ymax7_n6p5.mat'), 'step_ref');
%   load(fullfile(step_root, 'many_steps_data_rand_ymax7.mat'), 'step_ref');
  yref = step_ref.yref;
  step_descr = 'many_steps_ymax7';
end

rw = 2e-07;
rng(1);
thenoise = timeseries(mvnrnd(0, rw, length(yref.Time))*0, yref.Time);


F_yudu = figure(60); clf
subplot(3,1,1)
hold on, grid on;
step_ref.plot(F_yudu, '-k', 'LineWidth', 0.5)

F_y = figure(61); clf
hold on, grid on
if max(abs(yref.Data)) > 0
  step_ref.plot(F_y);
  step_ref.plot_settle_boundary(F_y, TOL, tol_mode);
end

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
du_max = du_max_orig/norm(plants.gdrift_inv, Inf);
du_max = 1000;


if md == 1 
  cmplx_rad = 0.9;
  [Q1, R0, S1, P_x] = build_control_constsigma(plants.sys_recyc, cmplx_rad);
  % gam_lin = 100;
  % gam_mpc = 100;
elseif md ==2
  can_cntrl = CanonCntrlParams_ns14();
  [Q1, R0, S1, P_x] = build_control(plants.sys_recyc, can_cntrl);
else
  error('expected md =(1|2|3|4) but recieved md = %d', md);
end

% ------------------------- Observer Gain ---------------------------------
can_obs_params = CanonObsParams_01();
can_obs_params.beta = 50;
[sys_obsDist, L_dist] = build_obs(plants.SYS, can_obs_params);

% 2). Design FeedForward gains.
[Nx, Nu] = SSTools.getNxNu(plants.sys_recyc);

gams = logspace(log10(0.1), log10(1000), 100);
min_gam_slf_cs = 7.5;
min_gam_mpc_cs = 0.001;
min_gam_slf_cz = 3.5;
min_gam_mpc_cz = 0.00001;
gams = unique([gams, 1.0, 3.0, 3.1, 3.8, min_gam_slf_cs, min_gam_mpc_cs,...
  min_gam_slf_cz, min_gam_mpc_cz]);

gams = [0.001, 0.01, 0.1, 1, 10, 25, 50.94, 75, 100, 200];
dist_traj = yref;

GM_s = gams*0;
PM_s = gams*0;
MagS_s = gams*0;
BW_s = gams*0;
TS_s = zeros(length(gams), 2);

g_lpf = zpk([], [0.95], 1, Ts);
g_lpf = g_lpf/dcgain(g_lpf);

for j=1:2
  if j ==1
    dist_scale= .25;
    g_gain = 1;
  else
    dist_scale = 0;
    g_gain = 1.1;
  end
  figure(10 + j); clf
  
  for k=1:length(gams)
    gam_lin = gams(k);
    
    dist_traj.Data = +yref.Data*dist_scale
    thenoise.Data = -yref.Data*0;
    dist_traj.Data = lsim(g_lpf, dist_traj.Data, dist_traj.Time);
    
    g_pert = zpk([], [], g_gain, Ts);
    % -----------------------------------------------------
    K_lqr = dlqr(plants.sys_recyc.a, plants.sys_recyc.b, Q1, R0+gam_lin, S1);
    
    % [Sens, Hyd, Hyr, Hyeta, Loop] = ss_loops_delta_dist(plants.SYS, plants.sys_recyc,...
    %   sys_obsDist, K_lqr, L_dist);
    
    [GM, PM, Sens_gain, ~, BW] = gmpm_vs_gam_recyc_obs(plants.SYS, plants.sys_recyc,...
      sys_obsDist, Q1, R0, S1, L_dist, gam_lin);
    
    GM_s(k) = GM;
    PM_s(k) = PM;
    MagS_s(k) = Sens_gain;
    BW_s(k) = BW;
    
%     F_clbode = figure(25);clf; hold on, grid on
%     Hbode_sens = bodeplot(Sens);
%     setoptions(Hbode_sens, 'FreqUnits', 'Hz')
%     legend('S: linear', 'S: mpc')
%     grid on;
%     
%     
%     [Gm_lin, Pm_lin] = margin(Loop);
%     fprintf('-------- MARGINS ------------------\n')
%     fprintf('Linear: GM = %.2f [], PM = %.2f [deg]\n', Gm_lin, Pm_lin)

    subplot(2,2,1)
    yyaxis left
    set(gca(), 'XTick', [0.01, 1, 100])
    semilogx(gams(1:k), GM_s(1:k))
    title('GM/PM')
    yyaxis right
    semilogx(gams(1:k), PM_s(1:k))
    hold on
    
    subplot(2,2,2)
    semilogx(gams(1:k), MagS_s(1:k))
    title('Mag S')
    set(gca(), 'XTick', [0.01, 1, 100])
    
    subplot(2,2,3)
    semilogx(gams(1:k), BW_s(1:k))
    title('bandwidth')
    set(gca(), 'XTick', [0.01, 1, 100])
    
    % -------------------------------------------------------------------
    % --------------------  FP Linear stuff -----------------------------
    
    sims_fxpl = SimAFM(plants.PLANT, K_lqr, Nx, sys_obsDist, L_dist, du_max,...
      false, 'thenoise', thenoise);
    sims_fxpl.gdrift = g_pert;
    [y_fxpl, U_full_fxpl, U_nom_fxpl, dU_fxpl, Xhat_fxpl] = sims_fxpl.sim(yref, dist_traj);
    name = sprintf('FXP lin Sim. (%s)', exp_id_str);
    
    fxpl_Opts = stepExpDuOpts('pstyle', '-r', 'TOL', TOL, 'step_ref', step_ref,...
      'controller', K_lqr, 'name',  name);
    sim_exp_fxpl = stepExpDu(y_fxpl, U_full_fxpl, dU_fxpl, fxpl_Opts);
    
    
    TS_s(k,j) = sum(sim_exp_fxpl.settle_time(TOL, tol_mode, 1));
    
    
    h2 = plot(sim_exp_fxpl, F_yudu, 'umode', 'both');
    legend([h2(1)])
    figure(F_y)
    h22 = sim_exp_fxpl.ploty(F_y);
    legend([h22]);
    
    figure(10+j)
    subplot(2,2,4)
    semilogx(gams(1:k), TS_s(1:k, j));
    title('total settle-time')
    hold on
    set(gca(), 'XTick', [0.01, 1, 100])
  end
end


% return
%
% if saveon
%   save(fullfile(save_root, [step_descr, '_linfxp_sim_', exp_id_str, '_',...
%     datestr(now, 'mm-dd-yyyy'), exp_num, '_.mat']), 'sim_exp_fxpl');
%        
%   save(fullfile(save_root, [step_descr, '_mpcfxp_sim_', exp_id_str, '_',...
%     datestr(now, 'mm-dd-yyyy'), exp_num, '_.mat']), 'sim_exp_fxpm');
% end
% 












