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
dataOut_path    = fullfile(PATHS.step_exp, outputDataName);
% labview reads desired trajectory here
refTrajPath     = fullfile(PATHS.step_exp, refTrajName);
% location of the vi which runs the experiment.

% folder into which to save results. Point
% process_settle_time_data.m here...
experiment_directory = ['many_steps_data_rand_', date, '_01'];
step_exp_root = fullfile(PATHS.exp, 'step-exps');
[status, message ] = mkdir(step_exp_root, experiment_directory);
save_root = fullfile(step_exp_root, experiment_directory);
%%
fprintf('\n\n\n\n')
TOL = 14/512; % max volts by pixels
% TOL = .01;
tol_mode = 'abs';
% which simulations to run


do_sim_mpcfxp = true;
do_sim_hyst = true;
do_inv_hyst = true;
do_drift = true;
do_invdrift = true;

plotstate = false;
plotpoles = false;
plotdriftpoles = false;
plotdriftbode = false;
saveon = true;

% We need four modes in this script:
% 1). Constant sigma, minimum-gamma 
% 2). Constant sigma, rob-optimal
% 3). Choose zeta, minimum-gamma
% 4). Choose zeta, rob-optimal
%
% For all cases, drift and hysteresis compensation remains
% unchanged. 
md = 2;

% !!!!!! These gamas should match the data in table II !!!!!!
if md == 1
  % 1). Constant sigma, minimum-gamma 
  gam_mpc = 2;
  exp_id_str = 'const-sig-min-gam';
elseif md == 2
  % 2). Constant sigma, rob-optimal
  gam_mpc = 129.2;
  exp_id_str = 'const-sig-rob-opt';
elseif md ==3
% 3). Choose zeta, minimum-gamma
  gam_mpc = 0.2;
  exp_id_str = 'choose-zet-min-gam';
elseif md ==4
% 4). Choose zeta, rob-optimal
  gam_mpc = 61.4;
  exp_id_str = 'choose-zet-rob-opt';
elseif md == 5
  gam_mpc = 2;
  exp_id_str = 'const-sig-same-sig';
end

% ------- Load Plants -----
[plants, frf_data] = CanonPlants.plants_ns12;
Ts  = plants.SYS.Ts;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                  Design reference "trajectory"                          %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Get a ref trajectory to track.
N  = 800;
r1 =1.37;
r2 = 0;
trajstyle =2;
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
  load(fullfile(step_root, 'many_steps_data_rand_ymax7.mat'), 'step_ref');
  yref = step_ref.yref;
  step_descr = 'many_steps_ymax7';
end
dist_traj = yref;
dist_traj.Data = dist_traj.Data*0 + .25*0;
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
  % legend([h1(1)])
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
% Design control/estimator gains. This gains are what we actually         %
% use in both the simulation and experimentant system.                    %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% -------------------- Constrained LQR Stuff ------------------------------

% Adjust the du_max to account for the gain of gdrift_inv.
du_max_orig = StageParams.du_max;
du_max = du_max_orig/norm(plants.gdrift_inv, Inf);
% du_max = 1000;


if md == 1 || md == 2 || md == 5
  cmplx_rad = 0.9;
  [Q1, R0, S1, P_x] = build_control_constsigma_ns12(plants.sys_recyc, cmplx_rad);
elseif md ==3 || md ==4 
  can_cntrl = CanonCntrlParams_ns14();
  % % % % can_cntrl = can_cntrl.aggressive_params();
  [Q1, R0, S1, P_x] = build_control(plants.sys_recyc, can_cntrl);
else
  error('expected md =(1|2|3|4) but recieved md = %d', md);
end

K_mpc = dlqr(plants.sys_recyc.a, plants.sys_recyc.b, Q1, R0+gam_mpc, S1)*1;

N_mpc = 22;

Qp = dare(plants.sys_recyc.a, plants.sys_recyc.b, Q1, R0+gam_mpc, S1);
nu = 1;
mpcProb = condensedMPCprob_OA(plants.sys_recyc, N_mpc, Q1, Qp, R0+gam_mpc, S1);
% mpcProb.lb = zeros(N_mpc,1)-du_max;
% mpcProb.ub = zeros(N_mpc,1)+du_max;

CON = CondenCon([], [], N_mpc);
CON.add_input_con('box', [-du_max, du_max]);
% GI = ss(plants.G_uz2powI);
% CON = CondenCon(GI, GI.b*0, N_mpc);
% CON.add_state_con('box', [0.1])
mpcProb.CON = CON;

Hmpc = mpcProb.H; Mmpc = mpcProb.M;
maxIter = 20;
fprintf('condition of H = %.1f\n', mpcProb.kappa);

fgm_fp = FGMprob_1(plants.sys_recyc, N_mpc, Q1, Qp, R0+gam_mpc, S1, du_max, maxIter);
I_H_mpc = fgm_fp.I_HL;
ML_x0  = fgm_fp.ML;
beta    = fgm_fp.beta;
fprintf('Mx0 needs n_int = %d\n', ceil(log2(double(max(abs(ML_x0(:))))))+1);
fprintf('IHL needs n_int = %d\n', ceil(log2(double(max(abs(I_H_mpc(:))))))+1);

% ------------------------- Observer Gain ---------------------------------
can_obs_params = CanonObsParams_01();
can_obs_params.beta = 50;
[sys_obsDist, L_dist] = build_obs(plants.SYS, can_obs_params);

[Sens_mpc, Hyd_mpc, Hyr_mpc, Hyeta_mpc, Loop_mpc] = ss_loops_delta_dist(plants.SYS,...
  plants.sys_recyc, sys_obsDist, K_mpc, L_dist);

F_clbode = figure(25);clf; hold on, grid on
Hbode_sens = bodeplot(Sens_mpc);
setoptions(Hbode_sens, 'FreqUnits', 'Hz')
legend('S: mpc')
grid on;

[Gm_mpc, Pm_mpc] = margin(Loop_mpc);
fprintf('-------- MARGINS ------------------\n')
fprintf('MPC: GM = %.2f [], PM = %.2f [deg]\n', Gm_mpc, Pm_mpc)

figure(102);
rlocus(Loop_mpc);
title('Kmpc')

figure(103)
nyquist(Loop_mpc)
title('Kmpc')
xlim([-2, .1])
ylim([-1, 1])

% 2). Design FeedForward gains.
[Nx, Nu] = SSTools.getNxNu(plants.sys_recyc);

% ------------------------ MPC floating point simulation  ---------------
% (Deleted, if you need it, paste it in from FXP_linear_MPC_...)

% -------------------------------------------------------------------
% -------------------- Setup Fixed stuff -----------------------------
A_obs_cl = sys_obsDist.a - L_dist*sys_obsDist.c;
fprintf('A_cl needs n_int = %d\n', ceil(log2(max(max(abs(A_obs_cl))))) + 1);
fprintf('L needs n_int = %d\n', ceil(log2(max(abs(L_dist)))) + 1);
fprintf('Nx needs n_int = %d\n', ceil(log2(max(abs(Nx)))) + 1);
fprintf('B needs n_int = %d\n', ceil(log2(max(abs(sys_obsDist.b))))+2);

nw = 32;
nf = 26;

du_max_fxp = fi(du_max, 1, 32, 26);
Nx_fxp = fi(Nx, 1, 32, 30);
L_fxp = fi(L_dist, 1, 32, 30);

sys_obs_fxp.a = fi(sys_obsDist.a -L_dist*sys_obsDist.c, 1, nw, nw-7);
sys_obs_fxp.b = fi(sys_obsDist.b, 1, nw, 29);
sys_obs_fxp.c = fi(sys_obsDist.c, 1, nw, 28);

% ----------------------------------------------------------------------- %
% --------------------- MPC, fgm fixed-point ---------------------------- %
if do_sim_mpcfxp
  nw_fgm = 32;
  nf_fgm = 28;
  
  fgm_fxp = FGMprob_fxp_1(plants.sys_recyc, N_mpc, Q1, Qp, R0+gam_mpc, S1, du_max,...
    maxIter, nw_fgm, nf_fgm);
  fgm_fxp.ML = fi(fgm_fp.ML, 1, 32, 26);
  fgm_fxp.x_nw = 32;
  fgm_fxp.x_nf = 27;
  fprintf('Mx0 needs n_int = %d\n', ceil(log2(max(abs(fgm_fp.ML(:)))))+1);
  fprintf('I_HL needs n_int = %d\n', ceil(log2(max(abs(fgm_fp.I_HL(:)))))+1);
  
  sims_fxpm = SimAFM(plants.PLANT, fgm_fxp, Nx_fxp, sys_obs_fxp, L_fxp, du_max_fxp,...
    true, 'nw', nw, 'nf', nf, 'thenoise', thenoise);
  
  if do_sim_hyst
    sims_fxpm.r = plants.hyst_sat.r;
    sims_fxpm.w = plants.hyst_sat.w;
    sims_fxpm.rp = fi(plants.hyst_sat.rp, 1, 16, 11);
    sims_fxpm.wp = fi(plants.hyst_sat.wp, 1, 16, 11);
    sims_fxpm.d = plants.hyst_sat.d;
    sims_fxpm.ws = plants.hyst_sat.ws;
    sims_fxpm.dp = fi(plants.hyst_sat.dp, 1, 16, 11);
    sims_fxpm.wsp = fi(plants.hyst_sat.wsp, 1, 16, 11);
    sims_fxpm.gdrift_inv = plants.gdrift_inv;
    sims_fxpm.gdrift = plants.gdrift;
  end
  
  [y_fxpm, U_full_fxpm, U_nom_fxpm, dU_fxpm, Xhat_fxpm, Xerr_fxpm] = sims_fxpm.sim(yref, dist_traj);
  name = sprintf('FXP MPC Sim. (%s)', exp_id_str);
  fxpm_Opts = stepExpDuOpts('pstyle', '--g', 'TOL', TOL, 'step_ref', step_ref,...
    'controller', fgm_fxp, 'name',  name);
  
  sim_exp_fxpm = stepExpDu(y_fxpm, U_full_fxpm, dU_fxpm, fxpm_Opts);
  
  h3 = plot(sim_exp_fxpm, F_yudu, 'umode', 'both');
  legend([h3(1)]);
  
  h32 = sim_exp_fxpm.ploty(F_y);
  legend([h32]);
  
  
  Ts_vec_fxpm = sim_exp_fxpm.settle_time(TOL, tol_mode, 1);
  fprintf('Total MPC FXP settle-time = %.3f [ms]\n', sum(Ts_vec_fxpm)*1000);
  
  figure(1000);
  hold on
  Ipow = lsim(plants.G_uz2powI, U_full_fxpm.Data, U_full_fxpm.Time);
  sim_exp_fxpm.Ipow = timeseries(Ipow, U_full_fxpm.Time);
  
  hmpc_Ipow = plot(U_full_fxpm.Time, Ipow, '--g');
  hmpc_Ipow.DisplayName = 'MPC Current';
  
  fprintf('\n-------------------------------------------------\n')
  [ts_mat, names] = pretty_print_ts_data(TOL, tol_mode, sim_exp_fxpm);

end
% return

if saveon
  save(fullfile(save_root, [step_descr, '_mpcfxp_sim_', exp_id_str, '_',...
    datestr(now, 'mm-dd-yyyy'), '02.mat']), 'sim_exp_fxpm');
end

sims_fxpm.sys_obs_fp = sys_obsDist;
sims_fxpm.sys_obs_fp.a = sys_obsDist.a - L_dist*sys_obsDist.c;


traj_path = 'Z:\mpc-journal\step-exps\traj_data.csvtraj_data.csv';

%%
%--------------------------------------------------------------------------
% -------------- MPC Experiment -------------------------------------------

% Build the u-reset.
if 1
  dry_run = false;
  reset_piezo('t1', 15, 't_final', 25, 'umax', 9, 'k1', 0.55,...
            'verbose', true, 'dry_run', dry_run)
end
%%
if ispc & do_sim_mpcfxp
  mpc_dat_path = 'Z:\mpc-journal\step-exps\MPCControls01.csv';
  sims_fxpm.write_control_data(mpc_dat_path, yref, traj_path)
  'pp'
end
%
SettleTicks = 20000;
Iters = length(yref.Data)-1;
% Iters = 200;
% create and pack data. Then save it.

[num, den] = tfdata(plants.gdrift_inv);
num = num{1};
den = den{1};

umax = 3.3;
ymax = max(yref.Data)*1.3
clear e;
clear vi;
%
% -----------------------RUN THE Experiment--------------------------------
vipath =['C:\Users\arnold\Documents\MATLAB\afm_mpc_journal',...
  '\labview\fixed-point-host\play_FXP_AFMss_MPC_distEst_singleAxisN22.vi'];
if 1
[e, vi] = setupVI(vipath, 'SettleTicks', SettleTicks, 'Iters', Iters,...
   'num', num, 'den', den, 'TF Order', 1*(length(den)-1),...
   'r_s', plants.hyst_sat.rp, 'w_s', plants.hyst_sat.wp, 'N_hyst', 1*length(plants.hyst_sat.rp),...
   'sat_ds', plants.hyst_sat.dp, 'sat_ws', plants.hyst_sat.wsp, 'N_sat', 1*length(plants.hyst_sat.dp),...
   'du_max', du_max,'dry_run', false,...
   'umax', umax, 'ymax', ymax, 'outputDataPath', dataOut_path,...
   'traj_path', traj_path, 'control_data_path', mpc_dat_path);
else
[e, vi] = setupVI(vipath, 'SettleTicks', SettleTicks, 'Iters', Iters,...
          'num', num, 'den', den, 'TF Order', (length(den)-1),...
          'r_s', plants.hyst.rp, 'w_s', plants.hyst.wp,'N_hyst', length(plants.hyst_sat.rp),...
            'umax', umax, 'ymax', ymax, 'outputDataPath', dataOut_path,...
            'dry_run', false,...
            'traj_path', traj_path, 'control_data_path', mpc_dat_path);
end          
vi.Run
%
ticks_bench = vi.GetControlValue('Loop Ticks (benchmark)');
fprintf('Actual loop ticks: %d\n', ticks_bench);
% Now, read in data, and save to structure, and plot.
AFMdata = csvread(dataOut_path);

t_exp = (0:size(AFMdata,1)-1)'*Ts;
y_exp = timeseries(AFMdata(:,1), t_exp);
u_exp = timeseries(AFMdata(:, 2), t_exp);
du_exp = timeseries(AFMdata(:,3), t_exp);
ufull_exp = timeseries(AFMdata(:,4), t_exp);

Ipow_exp = timeseries(AFMdata(:,5), t_exp);
xhat_exp = timeseries(AFMdata(:,6:end), t_exp);
yy = xhat_exp.Data*sys_obsDist.c';

expOpts = stepExpDuOpts('TOL', TOL, 'step_ref', step_ref, 'controller', fgm_fxp,...
  'pstyle', '--m', 'name', sprintf('AFM Stage (MPC, %s), ', exp_id_str));

afm_exp_mpc = stepExpDu(y_exp, ufull_exp, du_exp, expOpts);
afm_exp_mpc.Ipow = Ipow_exp;

try
Ts_vec_afm_mpc = afm_exp_mpc.settle_time(TOL, tol_mode, 1);
fprintf('Total AFM mpc FXP settle-time = %.3f [ms]\n', sum(Ts_vec_afm_mpc)*1000);
end

H_mpcexp = plot(afm_exp_mpc, F_yudu, 'umode', 'both');
try; legend([h1(1), h2(1), h3(1),  H_mpcexp(1) ]); end
subplot(3,1,1)
plot(y_exp.Time, yy, ':k')

H_mpcexp2 = afm_exp_mpc.ploty(F_y);

try; legend([h12, h22, h32, H_mpcexp2]); end

figure(1000); %clf
plot(Ipow_exp.Time, (Ipow_exp.Data/15)*1000)
ylabel('current [mA]')
grid on
title('Current')

try; [~, F_state] = plotState(xhat_exp, F_state); end
fprintf('Max of experimental Xhat = %.2f\n', max(abs(xhat_exp.data(:))));
[ts_mat, names] = pretty_print_ts_data(TOL, tol_mode,...
  sim_exp_fxpm, afm_exp_mpc);


% save('many_steps_data/many_steps_rand_fxpmpc_invHystDrift_constsig_7-9-2018.mat', 'y_exp', 'u_exp',...
%   'du_exp', 'ufull_exp', 'Ipow_exp', 'yref', 'y_fxpm')
save(fullfile(save_root, [step_descr, '_mpc_EXP_', exp_id_str, '_',...
  datestr(now, 'mm-dd-yyyy'), '.mat']), 'afm_exp_mpc');















