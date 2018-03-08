% I want to investigate methods to increase the condition number of the
% Hessian

% Build up the correct model from what is saved from sysID. Ie, put the
% thing in a 
clear
clc
close all
addpath('functions')
addpath('models')

volts2mu = 1;
TOL = 0.01;

fname_lin = 'data/max_sp_data_lin_TOY.mat';
fname_mpc = 'data/max_sp_data_mpc_TOY.mat';
fname_clqr = 'data/clqr_ref_data_TOY.mat';
fname_timeopt = 'data/timeopt_ref_data_TOY.mat';


% ----------------------------------------------------------------------- %
% ----------------------- Create a toy model ---------------------------- %
Ts = 40e-6;

w1 = 650*2*pi;
z1 = 0.02;
w2r = 500*2*pi;
s1 = -w1*z1 + 1j*w1*sqrt(1-z1^2);
s2 = conj(s1);
s3 = -w2r;
% real poles
ps = [s1, s2, s3];
% discrete poles
pz = exp(ps*Ts);
G = zpk([], pz, 1, Ts);
sys = ss(G)/dcgain(G);

[Gfrf, ws] = freqresp(sys, logspace(0, 3, 100));
Gfrf = Gfrf(:);
F1 = figure(1);
frfBode(sys, ws*2*pi, F1, '-r', 'Hz')
plotPZ_freqs(sys, F1)


sys_nodelay = sys;
PLANT = sys;
Ns = length(sys.b);
Nd = 0;

xss = SSTools.getNxNu(PLANT);
dcgain_sys = 1/(PLANT.c*xss);
x0 = xss*0;


% 3). LQR generation gain.        
% -------------------------------------------------------------------------
% -------------------- Constrained LQR Stuff ------------------------------



zeta_x = [.9];
gams_x = [2.5];
rhos_x = 2.5;
 

P_x    = getCharDes(sys_nodelay, gams_x, [], zeta_x, rhos_x, .25);
K_temp = place(sys_nodelay.a, sys_nodelay.b, P_x);
[Q0, R1, K_lqr] = inverseLQR(sys_nodelay, K_temp);

R1 = 1000;

sys_recyc = SSTools.deltaUkSys(sys);
Ns_mpc = size(sys_recyc.B, 1);
[Nx, Nu] = SSTools.getNxNu(sys_recyc);

Q1 = blkdiag(Q0, 0);
K_lqr = dlqr(sys_recyc.a, sys_recyc.b, Q1, R1);
Qp = dare(sys_recyc.a, sys_recyc.b, Q1, R1); 


% ------------------ Linear + delU Saturation Case ---------------------- %
% Iterate over a bunch of gammas and try to find the maximum setpoint for
% each one. This code is pretty niave. We start at a very low setpoint and
% slowly increase the setpoint until the settling time is reported as NaN.
% A bisection search would probably be faster. This is complicated though
% by the fact that instability seems to occur in the "middle": if the
% setpoint is large enough, we dont have the stability problem, which is
% weird.
% ref_s = linspace(0.01, 15, 200);

R = 5;
mpc_on=0;

%%

% *Optimal, open-loop* trajectory generation over all setpoints. This
% provides the basiline for comparison of the other methods. 
% ref_s = [1.5];
du_max   = 0.5;
gamma = 1;
gam_s = linspace(gamma, 100, 10); % original
% gam_s = [1, 100, 1000, 2500, 5000, 10000];
ref_s = 0.1:0.25:10;


% N_mpc_s = [4, 8, 12, 16, 20]; % original 
N_mpc_s = [12];
N_traj =800;
trun = Ts*N_traj;

logfile = sprintf('logs/log-lin-mpc-parfor_TOY_%s.log', date);
LG = EchoFile(logfile);
logger = @LG.echo_file;
logger = @fprintf;
ProgBar = @(max_iter, varargin)ProgressBarFile(max_iter, varargin{:}, 'logger', logger);
%
% warning('ON', 'MATLAB:mir_warning_maybe_uninitialized_temporary');
step_params_lin = StepParamsLin(sys_recyc, ref_s, du_max,Q1, gam_s, PLANT, trun);
step_data_lin = StepDataLin(step_params_lin, 'savedata', true,...
    'file', fname_lin', 'logger', logger,...
    'Progbar', ProgBar);

step_params_mpc = StepParamsMPC(sys_recyc, ref_s, du_max,Q1, gam_s, PLANT,...
                    trun, N_mpc_s,'condensed');
step_data_mpc = StepDataMPC(step_params_mpc, 'savedata', true,...
    'file', fname_mpc', 'logger', logger,...
    'Progbar', ProgBar);


% step_params_timeopt = StepParamsTimeOpt(sys, ref_s, du_max, sys_nodelay, 10);
step_params_clqr  = StepParamsCLQR(sys_recyc, ref_s, du_max,Q1, gamma, PLANT, N_traj, 'condensed');
step_data_clqr = StepDataCLQR(step_params_clqr, 'savedata', true,...
    'file', fname_clqr);

step_params_timeopt = StepParamsTimeOpt(sys_recyc, ref_s, du_max, sys_nodelay, Nd);
step_data_timeopt = StepDataTimeOpt(step_params_timeopt, 'savedata', true,...
    'file', fname_timeopt);

%
% ======================================================================= %
%                                                                         %
%                        BUILD TRAJECTORIES                               %
%                                                                         %
% ======================================================================= %

% ------------- Generate/Load CLQR max Trajectories --------------------- %
step_data_timeopt = step_data_timeopt.build_timeopt_trajs('force', 0,...
                    'verbose', 3, 'max_iter', 50, 'do_eject', 0);

%
% ------------- Generate/Load CLQR max Trajectories --------------------- %
tic
try
    
    step_data_clqr = build_clqr_trajs(step_data_clqr, 'force', 0, 'verbose', 3);
    logger('Finished building clqr_data. Total time = %.2f\n', toc);
catch ME
    errMsg = getReport(ME, 'extended', 'hyperlinks', 'off');
    logger('Failed to build clqr_data: \n%s', errMsg);
    fprintf('Failed to build clqr_data: \n%s', errMsg);
end


ref_s = step_data_clqr.params.ref_s;
gam_s = step_data_clqr.params.gam_s;
ref_s_to = step_data_timeopt.params.ref_s;

%%
F200=figure(200);clf, hold on;
colrs =  get(gca, 'colororder');
ax = gca;

hcl = step_data_clqr.plot_ref_vs_settle(ax,[], 'LineWidth', 1);
hto = step_data_timeopt.plot_ref_vs_settle(ax, 'LineWidth', 2);
legend([hcl, hto]);
%%
f10 = figure(10); clf
kk = 25
ax1= subplot(2,1,1); hold on, grid on;
ax2 = subplot(2,1,2); hold on, grid on;
step_data_clqr.plot_single_traj(kk, ax1, ax2)
step_data_timeopt.plot_single_traj(kk, ax1, ax2)


linkaxes([ax1, ax2], 'x')
xlim([0, 0.005])

f10.CurrentAxes = ax1;
ref_k = ref_s(kk)
xlm = xlim;
plot(xlm, [ref_k*1.01, ref_k*1.01], ':k')
plot(xlm, [ref_k*0.99, ref_k*0.99], ':k')


sys_sim = SSTools.deltaUkSys(step_data_timeopt.params.sys_nodelay);
[Nx_sim, Nu_sim] = SSTools.getNxNu(sys_sim);
x0_sim = Nx_sim*0;

xf = Nx_sim*ref_k;
k0 = 37;
toBisect = TimeOptBisect(sys_sim, du_max);
toBisect.max_iter = 50;
toBisect.TOL = 1e-3;
[X, U, status]=toBisect.time_opt_bisect(x0_sim, xf, 'k0', k0);


u = [U.Data; ones(2*k0,1)*ref_k*Nu_sim]; % zero for deltaUk
t = [0:1:length(u)-1]'*sys_sim.Ts;

[y, t, x] = lsim(sys_sim, u, t, x0_sim);


plot(t, y, '--')
f10.CurrentAxes = ax2
plot(t, u, '--')

%%
% ------------- Generate LIN max setpoints ------------------------------ %
threshold = 125; % percent
Judge = MaxSpJudgeCLQR(step_data_clqr, threshold);

try
    
    step_data_lin.max_ref_judge = Judge;
    step_data_lin = step_data_lin.build_max_setpoints('force', 0, 'verbose', 3);
    
    logger('Finished building max setpoints, linear. Total time = %.2f\n\n', toc);
catch ME
    errMsg = getReport(ME, 'extended', 'hyperlinks', 'off');
     logger('Failed to build max setpoints, linear: \n\n%s', errMsg);
     fprintf('Failed to build max setpoints, linear: \n\n%s', errMsg);
 end
%%
% ---------------------- Generate  mpc max setpoints -------------------- %

tic
try
    step_data_mpc.max_ref_judge = Judge;
    step_data_mpc = step_data_mpc.build_max_setpoints('force', 0, 'verbose', 1);
   logger('Finished building max setpoints, mpc. Total time = %.2f\n\n', toc);
catch ME
    errMsg = getReport(ME,  'extended','hyperlinks', 'off');
    logger('Failed to build max setpoints, mpc: %s\n\n', errMsg); 
    fprintf('Failed to build max setpoints, mpc: %s\n\n', errMsg); 
end

%%
% expose variables for plotting.
clear StepDataCLQR
clear StepDataTimeOpt
ref_s = step_data_clqr.params.ref_s;
gam_s = step_data_clqr.params.gam_s;
ref_s_to = step_data_timeopt.params.ref_s;


F200=figure(200);clf, hold on;
colrs =  get(gca, 'colororder');
ax = gca;

hcl = step_data_clqr.plot_ref_vs_settle(ax,[], 'LineWidth', 1);
hto = step_data_timeopt.plot_ref_vs_settle(ax, 'LineWidth', 2);
legend([hcl, hto]);
saveon = 0;

if saveon
    saveas(F200, 'figures/clqrTimeOpt_sp_vs_ts_CCTA.svg')
end
%%
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
legend(hands_mxsp)

xlabel('$\gamma$', 'interpreter', 'latex', 'FontSize', 14);
ylabel('Max Ref', 'interpreter', 'latex', 'FontSize', 14);


%%
plot(ref_s, step_data_clqr.results.settle_times_opt_cell{1})

%%
% ======================================================================= %
%                                                                         %
% ------------------Settle-By-Maximum-Acheivable-Reference -------------- %
%                                                                         %
% ======================================================================= %

% ---------------------------- For MPC ---------------------------------- %
rmax_s = [1, 2.5, 5.0, 10];
for jj = 1:length(N_mpc_s)
    N_mpc = N_mpc_s(jj);
    F=figure(300 + jj);clf; hold on;
    colrs =  get(gca, 'colororder');
    ax = gca;

    
    hto = step_data_timeopt.plot_ref_vs_settle(ax, 'LineWidth', 2,...
                            'Color', colrs(1, :));
    hcl = step_data_clqr.plot_ref_vs_settle(ax,[], 'LineWidth', 2,...
                         'Color', colrs(2, :));
    hands = [hcl; hto];

    
    for kk = 1:length(rmax_s)
        step_data_mpc = step_data_mpc.ts_by_ref_max(rmax_s(kk), jj);
        figure(F);
        [~, h] = step_data_mpc.plot_ts_by_ref_max_judge(rmax_s(kk),...
                       jj, ax, 'LineWidth', 1.25, 'Color', colrs(kk + 2,:));
        hands = [hands; h];
    end
    title(sprintf('N = %.0f', N_mpc))
    legend(hands);
end

% ------------------------- For LINEAR ---------------------------------- %
clc
F=figure(300);clf; hold on;
colrs =  get(gca, 'colororder');
ax = gca;
hto = step_data_timeopt.plot_ref_vs_settle(ax, 'LineWidth', 2,...
                        'Color', colrs(1, :));
hcl = step_data_clqr.plot_ref_vs_settle(ax,[], 'LineWidth', 2,...
                     'Color', colrs(2, :));
hands = [hcl; hto];
for kk = 1:length(rmax_s)
    step_data_lin = step_data_lin.ts_by_ref_max(rmax_s(kk), 1);
    figure(F)
    [~, h] = step_data_lin.plot_ts_by_ref_max_judge(rmax_s(kk),...
                   [], ax, 'LineWidth', 1.25, 'Color', colrs(kk + 2,:));
    hands = [hands; h];
end
legend(hands);
title('Linear')


%
% ======================================================================= %
%                                                                         %
%               Settle Time Percent Increase over time-Optimal            %
%                       by-Maximum-Acheivable-Reference                   %
%                                                                         %
% ======================================================================= %
clear StepDataCLQR
ts_timeopt = step_data_timeopt.results.settle_times_opt_cell{1};

% -----------------------CLQR vs Time-Optimal---------------------------- %
F = figure(399);
ax = gca();
step_data_clqr.plot_ts_perc_increase_by_rmax(1, ts_timeopt, ax);
title('CLQR')


%%% plot(ref_s(kref), (ts_mpc/ts_clqr)*100, 'xk')
% ------------------------- For LINEAR ---------------------------------- %
clc
clear StepDataQuad
F=figure(400);clf; hold on;
ax = gca();

hands = [];
for kk = 1:length(rmax_s)
    figure(F)
    step_data_lin = step_data_lin.ts_by_ref_max(rmax_s(kk), 1);
    [~, h] = step_data_lin.plot_ts_perc_increase_by_rmax(rmax_s(kk),...
                   1, ts_timeopt, ax, 'LineWidth', 1.25, 'Color', colrs(kk, :));
    hands = [hands; h];

end
title('Linear (Percent Increase)')
legend(hands)

for jj = 1:length(N_mpc_s)
    N_mpc = N_mpc_s(jj);
    F=figure(400 + jj);clf; hold on;
    colrs =  get(gca, 'colororder');
    ax = gca;

    hands = [];
    
    for kk = 1:length(rmax_s)
        step_data_mpc = step_data_mpc.ts_by_ref_max(rmax_s(kk), jj);
        figure(F);
        [~, h] = step_data_mpc.plot_ts_perc_increase_by_rmax(rmax_s(kk),...
                       jj, ts_timeopt, ax, 'LineWidth', 1.25, 'Color', colrs(kk + 2,:));
        hands = [hands; h];
    end
    title(sprintf('N = %.0f (Percent Increase)', N_mpc))
    legend(hands);
    grid on; zoom on;
end
%
clc
clear StepDataCLQR
clear StepDataQuad
for kk = 1:length(rmax_s)
    
   F = figure(500 + kk); clf; hold on;
   ax = gca();
   hands = [];
   
   % CLQR: this guy doesn't change, no need to regenerate.
   [~, h] = step_data_clqr.plot_ts_perc_increase_by_rmax(1, ts_timeopt, ax, 'Color', colrs(1, :));
   h.DisplayName = 'CLQR';
   hands = [hands;h];
   
   
   % -- Linear
   step_data_lin = step_data_lin.ts_by_ref_max(rmax_s(kk), 1);
   [~, h] = step_data_lin.plot_ts_perc_increase_by_rmax(rmax_s(kk),...
                   1, ts_timeopt, ax, 'LineWidth', 1.5, 'Color', colrs(2, :));
   h.DisplayName = sprintf('Linear ($\\gamma = %.0f$)', step_data_lin.ts_by_rmax_results{1}.gamma);
   hands = [hands;h];
   
   step_data_mpc.ts_by_rmax_results = {};
   for jj = 1:length(N_mpc_s)
        step_data_mpc = step_data_mpc.ts_by_ref_max(rmax_s(kk), jj);
        figure(F);
        [~, h] = step_data_mpc.plot_ts_perc_increase_by_rmax(rmax_s(kk),...
                       jj, ts_timeopt, ax, 'LineWidth', 1.5, 'Color', colrs(jj + 2,:));
        h.DisplayName = sprintf('MPC ($N=%.0f$ $\\gamma = %.0f$', N_mpc_s(jj),...
                        step_data_mpc.ts_by_rmax_results{jj}.gamma);
        hands = [hands; h];
    end
   
   legend(hands)
   title(sprintf('r-max = %.2f', rmax_s(kk)));
    
end
%%
tilefigs(1)
%%
          
group0 =  [figure(300), figure(301), figure(302), figure(303), figure(304),...
          figure(305)];
tilefigs([1], group0);
%%
group1 = [figure(400), figure(401), figure(402), figure(403), figure(404),...
          figure(405)];
tilefigs(2, group1)      

%%
group2 = [figure(501), figure(502),figure(503),figure(504)];      
tilefigs(3, group2);




