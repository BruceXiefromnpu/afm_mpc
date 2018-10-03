% I want to investigate methods to increase the condition number of the
% Hessian

% Build up the correct model from what is saved from sysID. Ie, put the
% thing in a
clear
clc
close all
addpath('functions')
addpath('models')
%

clc
fname_lin = fullfile(PATHS.sim_data(), 'max_sp_data_lin_dumax_p198.mat');
fname_mpc = fullfile(PATHS.sim_data(),'max_sp_data_mpc_dumax_p198.mat');
fname_clqr = fullfile(PATHS.sim_data(),'clqr_ref_data__dumax_p198.mat');
fname_timeopt = fullfile(PATHS.sim_data(),'timeopt_ref_data_dumax_p198.mat');

matpath           = getMatPath();
dataroot          = fullfile(matpath, 'afm_mpc_journal', 'modelFitting', 'pow_amp');

% volts2mu = 1;
TOL = 0.01;
umax = 5;

% 3). LQR generation gain.
% -------------------------------------------------------------------------
% -------------------- Constrained LQR Stuff ------------------------------
% du_max   = StageParams.du_max_vib;
% plants = CanonPlants.plants_with_drift_inv(false);
% sys_recyc = plants.sys_recyc;
% Ts = sys_recyc.Ts;
% sys_recyc_nodelay = plants.sys_recyc_nodelay;
% 

plants = CanonPlants.plants_ns14(9,2);
du_max_orig = StageParams.du_max;
du_max = du_max_orig/norm(plants.gdrift_inv, Inf);


sys_recyc = plants.sys_recyc;
Ts = sys_recyc.Ts;
sys_recyc_nodelay = plants.sys_recyc_nodelay;

can_cntrl = CanonCntrlParams_ns14()
[Q1, R0, S1] = build_control(sys_recyc, can_cntrl);

gam_lin = 3;
gam_mpc = 0.00001;
R1 = R0 + gam_mpc;
K_lqr = dlqr(sys_recyc.a, sys_recyc.b, Q1, R0+gam_lin, S1);
sys_cl = SSTools.close_loop(sys_recyc, K_lqr);

figure(2); clf
[y,t,x] = step(sys_cl, 800*Ts);
subplot(2,1,1)
plot(t,y)
grid on

subplot(2,1,2)
plot(t, (x*K_lqr'))
grid on


% ------------------ Linear + delU Saturation Case ---------------------- %
% Iterate over a bunch of gammas and try to find the maximum setpoint for
% each one. This code is pretty niave. We start at a very low setpoint and
% slowly increase the setpoint until the settling time is reported as NaN.
% A bisection search would probably be faster. This is complicated though
% by the fact that instability seems to occur in the "middle": if the
% setpoint is large enough, we dont have the stability problem, which is
% weird.


% *Optimal, open-loop* trajectory generation over all setpoints. This
% provides the basiline for comparison of the other methods.
% ref_s = [1.5];

R0 = 1;
gamma = 0.00001;
gam_s = linspace(gamma, 6, 30)
% gam_s = logspace(log10(gamma), log10(20), 30)
%
% gam_s = [1, 100, 1000, 2500, 5000, 10000];
ref_s = 0.1:0.5:15.5;


N_mpc_s = [8, 12, 22]; % original
% N_mpc_s = [12, 18, 24];
N_traj =800;
trun = Ts*N_traj;

clear StepDataMPC

logfile = sprintf('logs/log-lin-mpc-parfor_judge_dynamic_%s.log', date);
LG = EchoFile(logfile);
logger = @LG.echo_file;
logger = @fprintf;
ProgBar = @(max_iter, varargin)ProgressBarFile(max_iter, varargin{:}, 'logger', logger);

% warning('ON', 'MATLAB:mir_warning_maybe_uninitialized_temporary');

step_params_lin = StepParamsLin(sys_recyc, ref_s, du_max,Q1, R0, gam_s, plants.PLANT, trun, 'S', S1);

figs(1) = figure(10);
figs(2) = figure(11);
step_params_lin.sim(2, R0+gam_s(2), figs);


step_data_lin = StepDataLin(step_params_lin, 'savedata', true,...
    'file', fname_lin', 'logger', logger,...
    'Progbar', ProgBar);

step_params_mpc = StepParamsMPC(sys_recyc, ref_s, du_max,Q1, R0, gam_s, plants.PLANT,...
                    trun, N_mpc_s,'condensed', 'S', S1);
step_data_mpc = StepDataMPC(step_params_mpc, 'savedata', true,...
    'file', fname_mpc', 'logger', logger,...
    'Progbar', ProgBar);


% step_params_timeopt = StepParamsTimeOpt(sys, ref_s, du_max, sys_nodelay, 10);
step_params_clqr  = StepParamsCLQR(sys_recyc, ref_s, du_max,Q1, R0, gamma,...
  plants.PLANT, N_traj, 'condensed', 'S', S1);
step_data_clqr = StepDataCLQR(step_params_clqr, 'savedata', true,...
    'file', fname_clqr);

step_params_timeopt = StepParamsTimeOpt(sys_recyc, ref_s, du_max, sys_recyc_nodelay, plants.Nd);
step_data_timeopt = StepDataTimeOpt(step_params_timeopt, 'savedata', true,...
    'file', fname_timeopt);



% ======================================================================= %
%                                                                         %
%                        BUILD TRAJECTORIES                               %
%                                                                         %
% ======================================================================= %
verbose = 0;
close all
% ------------- Generate/Load Time-Opt max Trajectories ----------------- %
step_data_timeopt = step_data_timeopt.build_timeopt_trajs('force', 0,...
                    'verbose', verbose, 'max_iter', 50, 'TOL', 1e-4, 'do_eject', false);

% ------------- Generate/Load CLQR max Trajectories --------------------- %
tic
try
    step_data_clqr = build_clqr_trajs(step_data_clqr, 'force', 0, 'verbose', verbose);
    logger('Finished building clqr_data. Total time = %.2f\n', toc);
catch ME
    errMsg_logger = getReport(ME, 'extended', 'hyperlinks', 'off');
    errMsg = getReport(ME, 'extended', 'hyperlinks', 'on');
    logger('Failed to build clqr_data: \n%s', errMsg_logger);
    fprintf('Failed to build clqr_data: \n%s', errMsg);
end

% ------------- Generate LIN max setpoints ------------------------------ %
threshold = 120; % percent
Judge = MaxSpJudgeCLQR(step_data_clqr, threshold);

try

    step_data_lin.max_ref_judge = Judge;
    step_data_lin = step_data_lin.build_max_setpoints('force', 0, 'verbose', verbose);

    logger('Finished building max setpoints, linear. Total time = %.2f\n\n', toc);
catch ME
    errMsg = getReport(ME, 'extended', 'hyperlinks', 'off');
     logger('Failed to build max setpoints, linear: \n\n%s', errMsg);
     fprintf('Failed to build max setpoints, linear: \n\n%s', errMsg);
 end
%
% ---------------------- Generate  mpc max setpoints -------------------- %

tic
try
  
    step_data_mpc.max_ref_judge = Judge;
    step_data_mpc = step_data_mpc.build_max_setpoints('force', 0, 'verbose', verbose);
   logger('Finished building max setpoints, mpc. Total time = %.2f\n\n', toc);
catch ME
    errMsg = getReport(ME,  'extended','hyperlinks', 'off');
    logger('Failed to build max setpoints, mpc: %s\n\n', errMsg);
    fprintf('Failed to build max setpoints, mpc: %s\n\n', errMsg);
end

%%
% expose variables for plotting.
rmpath(genpath('~/matlab/solvers/cvx'));
ref_s = step_data_clqr.params.ref_s;
gam_s = step_data_clqr.params.gam_s;
ref_s_to = step_data_timeopt.params.ref_s;


F200=figure(200);clf; hold on;
colrs =  get(gca, 'colororder');
ax = gca;

hcl = step_data_clqr.plot_ref_vs_settle(ax,[], 'LineWidth', 2);
hto = step_data_timeopt.plot_ref_vs_settle(ax, 'LineWidth', 2);
legend([hcl, hto]);
grid on
saveon = 0;
%
if saveon
    saveas(F200, 'latex/figures/clqrTimeOpt_sp_vs_ts_CCTA_dumaxp6.svg')
end

% ----------------- Plot maximum reference vs gamma -----------------------
F10 = figure(10); clf; hold on;
colrs =  get(gca, 'colororder');
hands_mxsp = []
hands_mxsp(1) = plot(step_data_lin.params.gam_s, step_data_lin.results{1}.max_setpoints,...
                     '-', 'Color', colrs(1, :),...
            'LineWidth', 2);
set(hands_mxsp(1), 'DisplayName', 'LQR Lin + Sat');        
for k = 1:length(step_data_mpc.params.N_mpc_s)

    hands_mxsp(k+1) = plot(step_data_mpc.params.gam_s, step_data_mpc.results{k}.max_setpoints,...
                       '-', 'Color', colrs(k+1, :), 'LineWidth', 2);

    stit = sprintf('MPC, N=%.0f', step_data_mpc.params.N_mpc_s(k));
    set(hands_mxsp(k+1), 'DisplayName', stit);
end
leg1 = legend(hands_mxsp);
set(leg1, 'Position', [0.6775 0.4829 0.2181 0.2213]);

xlabel('$\gamma$', 'interpreter', 'latex', 'FontSize', 14);
ylabel('Max Ref', 'interpreter', 'latex', 'FontSize', 14);
grid on
if saveon
    saveas(F10, 'latex/figures/maxref_vs_gamma_dumax_Ipow_lowgain.svg')
end

%%
% ======================================================================= %
%                                                                         %
% ------------------Settle-By-Maximum-Acheivable-Reference -------------- %
%                                                                         %
% ======================================================================= %

% ---------------------------- For MPC ---------------------------------- %
F=figure(301);clf; hold on;
colrs =  get(gca, 'colororder');
ax = gca;
rmax_s = [14];
% rmax_s = [1, 2.5, 5.0, 14];
hto = step_data_timeopt.plot_ref_vs_settle(ax, 'LineWidth', 2,...
  'Color', colrs(1, :));
hcl = step_data_clqr.plot_ref_vs_settle(ax,[], 'LineWidth', 2,...
  'Color', colrs(2, :));
hands = [hcl; hto];

step_data_mpc = step_data_mpc.ts_by_ref_max(rmax_s, 1);
for jj = 1:length(N_mpc_s)
    N_mpc = N_mpc_s(jj);
    step_data_mpc = step_data_mpc.ts_by_ref_max(rmax_s, jj);
    gam_ = step_data_mpc.ts_by_rmax_results{jj}.gamma;
    step_data_mpc = step_data_mpc.ts_by_ref_max(rmax_s, jj);
    figure(F);
    [~, h] = step_data_mpc.plot_ts_by_ref_max_judge(rmax_s,...
      jj, ax, 'LineWidth', 1.25, 'Color', colrs(jj + 2,:));
    h.DisplayName = sprintf('MPC, $\\gamma=%.3g$, N=%d', gam_, N_mpc);
    hands = [hands; h];

    title(sprintf('N = %.0f', N_mpc))
    
end
legend(hands);

% ------------------------- For LINEAR ---------------------------------- %

step_data_lin = step_data_lin.ts_by_ref_max(rmax_s, 1);
gam_ = step_data_lin.ts_by_rmax_results{1}.gamma;
[~, h] = step_data_lin.plot_ts_by_ref_max_judge(rmax_s,...
  [], ax, 'LineWidth', 1.25, 'Color', colrs(jj + 3,:));
h.DisplayName = sprintf('Linear, $\\gamma=%.3g$, ', gam_);
hands = [hands; h];

legend(hands);
title('Linear')



ts_timeopt = step_data_timeopt.results.settle_times_opt_cell{1};

F = figure(501); clf; hold on;
ax = gca();
hands = [];

% CLQR: this guy doesn't change, no need to regenerate.
[~, h] = step_data_clqr.plot_ts_perc_increase_by_rmax(1, ts_timeopt, ax, 'Color', colrs(1, :));
h.DisplayName = 'CLQR';
hands = [hands;h];


% -- Linear
step_data_lin = step_data_lin.ts_by_ref_max(rmax_s, 1);
[~, h] = step_data_lin.plot_ts_perc_increase_by_rmax(rmax_s,...
  1, ts_timeopt, ax, 'LineWidth', 1.5, 'Color', colrs(2, :));
h.DisplayName = sprintf('Linear ($\\gamma = %.1f$)', step_data_lin.ts_by_rmax_results{1}.gamma);
hands = [hands;h];

step_data_mpc.ts_by_rmax_results = {};
for jj = 1:length(N_mpc_s)
  step_data_mpc = step_data_mpc.ts_by_ref_max(rmax_s, jj);
  figure(F);
  [~, h] = step_data_mpc.plot_ts_perc_increase_by_rmax(rmax_s,...
    jj, ts_timeopt, ax, 'LineWidth', 1.5, 'Color', colrs(jj + 2,:));
  h.DisplayName = sprintf('MPC ($N=%.0f$ $\\gamma = %.1f$', N_mpc_s(jj),...
    step_data_mpc.ts_by_rmax_results{jj}.gamma);
  hands = [hands; h];
end
grid on
legend(hands)
title(sprintf('r-max = %.2f', rmax_s));


%%
idx = find(step_data_lin.results{1}.max_setpoints >= rmax_s(end), 1, 'first');

ts_sums = [];
for k=idx:length(gam_s)
  
  t_settle_s = step_data_lin.results{1}.data{k}.t_settle_s;
  
  ts_sums(k-idx+1) = sum(t_settle_s.*ref_s);
end

gams_ = gam_s(idx:end);

figure(10); clf;
plot(gams_, ts_sums)
hold on; grid on;
%

idx = find(step_data_mpc.results{1}.max_setpoints >= rmax_s(end), 1, 'first');

ts_sums = [];
for k=idx:length(gam_s)
  
  t_settle_s = step_data_mpc.results{1}.data{k}.t_settle_s;
  
  ts_sums(k-idx+1) = sum(t_settle_s);
end

gams_ = gam_s(idx:end);

figure(10)
plot(gams_, ts_sums)
% cz_mg_mpc22 = step_data_mpc.ts_by_rmax_results{3}.gamma;
% cz_mg_slf   = step_data_lin.ts_by_rmax_results{1}.gamma;



% % gam = R1;
% % N = 20;
% % Qp = dare(sys_recyc.a, sys_recyc.b, Q1, gam, S1);
% % MPP1 = condensedMPCprob_OA(sys_recyc, N, Q1, Qp, gam, S1);
% % fprintf('kappa_fast = %.2f\n', MPP1.kappa)
%%
tilefigs(1)

figure(10)
xlm = xlim;
h1 = plot(xlm, [rmax_s', rmax_s'], ':k', 'LineWidth', 2)
if saveon
    saveas(F10, 'latex/figures/maxref_vs_gamma_dumax_Ipow_lowgain.svg')
end


%%

figure(504)
grid on
% xlim([0, 10])
saveas(gcf, fullfile(PATHS.jfig, 'perc_increase_Ipow_lowgain_rmax14.svg'))

% 
% figure(503)
% grid on
% % xlim([0, 5])
% saveas(gcf, fullfile(PATHS.jfig, 'perc_increase_Ipow_lowgain_rmax5.svg'))
% 
% figure(502)
% grid on
% % xlim([0, 2.5])
% saveas(gcf, fullfile(PATHS.jfig, 'perc_increase_Ipow_lowgain_rmax2p5.svg'))












