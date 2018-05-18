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

% folder into which to save results. Point
% process_settle_time_data.m here...
save_root = fullfile(PATHS.exp, 'experiments', 'many_steps_data');

% ---------- Load Parametric Models  -----------
% load(fullfile(PATHS.sysid, 'hysteresis/steps_hyst_model.mat'));
% load('C:\Users\arnold\Documents\MATLAB\afm_mpc_journal\modelFitting\hysteresis\steps_hyst_model_withsat.mat')
% load(fullfile(PATHS.sysid, 'FRF_data_current_stage2.mat'))
% r = r;
% w = theta_hyst;

umax = 5;

TOL = .01;
saveon = true;
%%
md = 2;
plants = CanonPlants.plants_with_drift_inv(false);

Ts  = plants.SYS.Ts;
if md == 2
  plants.gdrift = zpk([], [], 1, Ts);
  plants.gdrift_inv = zpk([], [], 1, Ts);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                  Design reference "trajectory"                          %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Get a ref trajectory to track.
N    = 800;
r1 = 6;
r2 = -6;
trajstyle =3;
if trajstyle == 1
  yref = CanonRefTraj.ref_traj_1(r1, N);
elseif trajstyle == 2
    yref = CanonRefTraj.ref_traj_2(r1, r2, N);
elseif trajstyle == 3
  yref = CanonRefTraj.ref_traj_load('many_steps.mat');
end
rw = 8.508757290909093e-07;
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
if md == 2
  du_max = du_max_orig/norm(plants.gdrift_inv);
else
  du_max = du_max_orig;
end


can_cntrl = CanonCntrlParams_01(plants.SYS);
[Q1, R0, S1] = build_control(plants.sys_recyc, can_cntrl);
gam_lin = 3;
gam_mpc = 0.5;
R1 = R0 + gam_mpc;

K_lqr = dlqr(plants.sys_recyc.a, plants.sys_recyc.b, Q1, R0+gam_lin, S1);
sys_cl = SSTools.close_loop(plants.sys_recyc, K_lqr);
N_mpc = 12;

Qp = dare(plants.sys_recyc.a, plants.sys_recyc.b, Q1, R1, S1);
nu = 1;
mpcProb = condensedMPCprob_OA(plants.sys_recyc, N_mpc, Q1, Qp, R1, S1);
% mpcProb.lb = zeros(N_mpc,1)-du_max;
% mpcProb.ub = zeros(N_mpc,1)+du_max;
CON = CondenCon([], [], N_mpc);
CON.add_input_con('box', [-du_max, du_max]);
mpcProb.CON = CON;

Hmpc = mpcProb.H; Mmpc = mpcProb.M;
maxIter = 20;
fprintf('condition of H = %.1f\n', mpcProb.kappa);

fgm_fp = FGMprob_1(plants.sys_recyc, N_mpc, Q1, Qp, R1, S1, du_max, maxIter);
I_H_mpc = fgm_fp.I_HL;
ML_x0  = fgm_fp.ML;
beta    = fgm_fp.beta;

% ------------------------- Observer Gain ---------------------------------
can_obs_params = CanonObsParams_01();
[sys_obsDist, L_dist] = build_obs(plants.SYS, can_obs_params);

if 1
  %%
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
  
  %saveas(F20, fullfile(PATHS.jfig, 'obs_cl.svg'))
    
end

% 2). Design FeedForward gains.
[Nx, Nu] = SSTools.getNxNu(plants.sys_recyc);


% for comparison, and for validation against the max_setpoint
% stuff, compute the time-optimal
sys_recyc_nod = SSTools.deltaUkSys(plants.sys_nodelay);
Nx_nod = SSTools.getNxNu(sys_recyc_nod);
tob = TimeOptBisect(sys_recyc_nod, du_max);
rmpath(genpath('~/matlab/solvers/cvx/lib'))
[xx_to, uu_to, stat] = tob.time_opt_bisect(sys_recyc_nod.b*0, Nx_nod*r1, 'k0', 50);
u = [uu_to.Data; zeros(100, 1)];
t = (0:length(u)-1)'*Ts;
yto = lsim(plants.sys_recyc, u, t);
ts_to = settle_time(t, yto, r1, TOL*r1);
figure, plot(t, yto)


sims_fpl = SimAFM(plants.PLANT, K_lqr, Nx, sys_obsDist, L_dist, du_max, false);
sims_fpm = SimAFM(plants.PLANT, mpcProb, Nx, sys_obsDist, L_dist, du_max, false);

if 1
%   sims_fpl.r = plants.hyst.r;
%   sims_fpl.w = plants.hyst.w;
%   sims_fpl.rp = plants.hyst.rp;
%   sims_fpl.wp = plants.hyst.wp;
  sims_fpl.gdrift_inv = plants.gdrift_inv;
  sims_fpl.gdrift = plants.gdrift;
end

[y_lin_fp_sim, U_full_fp_sim, U_nom_fp_sim, dU_fp_sim, Xhat_fp] = sims_fpl.sim(yref);
ts_lfp = settle_time(y_lin_fp_sim.Time(1:800), y_lin_fp_sim.Data(1:800), r1, TOL*r1);
fprintf('linear fp settle-time = %.3f [ms]\n', ts_lfp*1000);
fprintf('perc increase over time-optimal: %.3f\n', (ts_lfp/ts_to)*100);
linOpts = stepExpOpts('pstyle', '-r', 'TOL', TOL, 'y_ref', r1,...
                      'controller', K_lqr, 'name',  'FP lin Sim.');

sim_exp = stepExpDu(y_lin_fp_sim, U_full_fp_sim, dU_fp_sim, linOpts);

F1 = figure(59); clf
h1 = sim_exp.plot(F1, 'umode', 'both');
% H1 = plot(sim_exp, F1);
subplot(3,1,1)
hold on, grid on;
plot(yref.time, yref.Data, '--k', 'LineWidth', .05);
legend([h1(1)])


figure(70); clf
du_full = diff(U_full_fp_sim.Data);
du_full(end+1) = du_full(end);
hold on
plot(U_full_fp_sim.Time, du_full, '--g')
plot(dU_fp_sim.Time, dU_fp_sim.Data, 'r')

xlm = xlim();

plot(xlm, [du_max_orig, du_max_orig], ':k')
plot(xlm, -[du_max_orig, du_max_orig], ':k')
legend('du (actual)', 'du (nominal)')
grid on
F61 = figure(61); clf
plotState(Xhat_fp, F61);

if 1
  [y_fpm, U_fpm, ~, dU_fpm, Xhat_fpm, Xerr_fpm] = sims_fpm.sim(yref);
  fpm_Opts = stepExpOpts('pstyle', '--m', 'TOL', TOL, 'y_ref', r1,...
                          'controller', K_lqr, 'name',  'FP MPC Sim. (QP OA)');

  sim_exp_fpm = stepExpDu(y_fpm, U_fpm, dU_fpm, fpm_Opts);
  h3 = plot(sim_exp_fpm, F1, 'umode', 'both');

  legend([h1(1), h3(1)]);
  ts_mfp = settle_time(y_fpm.Time(1:800), y_fpm.Data(1:800), r1, TOL*r1);
  fprintf('mpc fp settle-time = %.3f [ms]\n', ts_mfp*1000);
  fprintf('perc increase over time-optimal: %.3f\n', (ts_mfp/ts_to)*100);

end

% -------------------- Setup Fixed stuff -----------------------------

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

% --------------------  Fixed Linear stuff -----------------------------

sims_fxpl = SimAFM(plants.PLANT, K_fxp, Nx_fxp, sys_obs_fxp, L_fxp, du_max_fxp,...
  true, 'nw', nw, 'nf', nf);
if 1
  % sims_fxpl.r = r;
  % sims_fxpl.w = w;
  % sims_fxpl.rp = rp;
  % sims_fxpl.wp = wp;
  sims_fxpl.gdrift_inv = plants.gdrift_inv;
  sims_fxpl.gdrift = plants.gdrift;
end


[y_fxpl, U_full_fxpl, U_nom_fxpl, dU_fxpl, Xhat_fxpl] = sims_fxpl.sim(yref);
fxpl_Opts = stepExpOpts('pstyle', '--k', 'TOL', TOL, 'y_ref', r1,...
                      'controller', K_lqr, 'name',  'FXP lin Sim.');
sim_exp_fxpl = stepExpDu(y_fxpl, U_full_fxpl, dU_fxpl, fxpl_Opts);

h2 = plot(sim_exp_fxpl, F1, 'umode', 'both');
legend([h1(1), h2(1), h3(1)])


[~, F61] = plotState(Xhat_fxpl, F61, [], [], '--');
fprintf('max of Xhat = %.2f\n', max(abs(Xhat_fxpl.Data(:))));
fprintf('max of M*Xhat = %.2f\n', max(max(abs(ML_x0*Xhat_fxpl.Data'))));



% --------------------- MPC, fgm fixed-point ------------------------
nw_fgm = 32;
nf_fgm = 28;

fgm_fxp = FGMprob_fxp_1(plants.sys_recyc, N_mpc, Q1, Qp, R1, S1, du_max,...
           maxIter, nw_fgm, nf_fgm);
fgm_fxp.x_nw = 32;
fgm_fxp.x_nf = 27;

sims_fxpm = SimAFM(plants.PLANT, fgm_fxp, Nx_fxp, sys_obs_fxp, L_fxp, du_max_fxp,...
                   true, 'nw', nw, 'nf', nf);

if 1
  [y_fxpm, U_fxpm, ~, dU_fxpm, Xhat_fxpm, Xerr_fxpm] = sims_fxpm.sim(yref);
  fxpm_Opts = stepExpOpts('pstyle', '--g', 'TOL', TOL, 'y_ref', r1,...
                          'controller', K_lqr, 'name',  'FXP MPC Simulation');

  sim_exp_fxpm = stepExpDu(y_fxpm, U_fxpm, dU_fxpm, fxpm_Opts);
  h3 = plot(sim_exp_fxpm, F1, 'umode', 'both');

  legend([h1(1), h2(1), h3(1)]);
  ts_mfxp = settle_time(y_fxpm.Time(1:800), y_fxpm.Data(1:800), r1, TOL*r1);
  fprintf('mpc fp settle-time = %.3f [ms]\n', ts_mfxp*1000);
end


if saveon
  save(fullfile(save_root, 'many_steps_linfp_sim.mat'),...
       'y_lin_fp_sim', 'U_full_fp_sim', 'U_nom_fp_sim', ...
       'dU_fp_sim')
  save(fullfile(save_root, 'many_steps_mpcfxp_sim.mat'),...
       'y_fxpm', 'U_fxpm', 'dU_fxpm');
end

return
sims_fxpm.sys_obs_fp = sys_obsDist;
sims_fxpm.sys_obs_fp.a = sys_obsDist.a - L_dist*sys_obsDist.c;

mpc_dat_path = 'C:\Users\arnold\Documents\MATLAB\afm_mpc_journal\data\MPCControls01.csv';
traj_path = 'C:\Users\arnold\Documents\MATLAB\afm_mpc_journal\data\traj_data.csv';
sims_fxpm.write_control_data(mpc_dat_path, yref, traj_path)


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
% Iters = 1200;
% Iters = min(Iters, length(yref)-1);
Iters = length(yref.Data);
% create and pack data. Then save it.

[num, den] = tfdata(gdrift_inv);
num = num{1};
den = den{1};

umax = 6.2;
ymax = max(yref.Data)*1.3
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
          'num', num, 'den', den, 'TF Order', (length(den)-1),...
          'r_s', rp, 'w-s', wp,'N_hyst', 7, 'dry_run', 0,...
            'umax', umax, 'ymax', ymax, 'outputDataPath', dataOut_path,...
            'traj_path', traj_path, 'control_data_path', mpc_dat_path);
          
vi.Run

ticks_bench = vi.GetControlValue('Loop Ticks (benchmark)');
fprintf('Actual loop ticks: %d\n', ticks_bench);
% Now, read in data, and save to structure, and plot.
AFMdata = csvread(dataOut_path);
% [y_exp, u_exp, I_exp, xhat_exp] = unpackExpData_nod(AFMdata, Ts);
t_exp = (0:size(AFMdata,1)-1)'*Ts;
y_exp = timeseries(AFMdata(:,1), t_exp);
u_exp = timeseries(AFMdata(:, 2), t_exp);
du_exp = timeseries(AFMdata(:,3), t_exp);
ufull_exp = timeseries(AFMdata(:,4), t_exp);
% xhat_exp = timeseries(AFMdata(:,5:end-1), t_exp);
Ipow_exp = timeseries(AFMdata(:,5), t_exp);
xhat_exp = timeseries(AFMdata(:,6:end), t_exp);

yy = xhat_exp.Data*sys_obsDist.c';
expOpts = stepExpOpts(linOpts, 'pstyle', '--m', 'name',  'AFM Stage (MPC)');

afm_exp = stepExpDu(y_exp, u_exp, du_exp, expOpts);
H2 = plot(afm_exp, F1, 'umode', 'both');
subplot(3,1,1)
plot(y_exp.Time, yy, ':k')

figure(1000); clf
plot(Ipow_exp.Time, (Ipow_exp.Data/15)*1000)
ylabel('current [mA]')
grid on
title('Current')

[~, F61] = plotState(xhat_exp, F61);
fprintf('Max of experimental Xhat = %.2f\n', max(abs(xhat_exp.data(:))));
%%
save('many_steps_data/many_steps_mpc_invHyst_invDrift2.mat', 'y_exp', 'u_exp', 'du_exp', 'wp', 'rp')






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
saveControlData(AllMatrix, 0, 0, Ns, Nd, r1, y_uKx, controlDataPath, refTrajPath);

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
% saveControlData(AllMatrix, 0, 0, Ns, Nd, r1, y_uKx, controlDataPath, refTrajPath);
%
clear e;
clear vi;
ns = length(K_lqr)-1;


dlmwrite(controlDataPath, AllMatrix, 'delimiter', ',', 'precision', 12)


Iters = 500
SettleTicks = 100;
% -----------------------RUN THE Experiment--------------------------------
vipath =['C:\Users\arnold\Documents\MATLAB\afm_mpc_journal\labview\',...
        'fixed-point-host\play_FXP_AFMss_LinearDistEst_singleAxis.vi'];
         
[e, vi] = setupVI(vipath, 'SettleTicks', SettleTicks, 'Iters', Iters,...
           'umax', 3, 'ymax', 5, 'du_max', du_max, 'Ns', ns, 'x_ref', r1,...
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






















