% I want to investigate methods to increase the condition number of the
% Hessian

% Build up the correct model from what is saved from sysID. Ie, put the
% thing in a 
clear
close all
addpath('functions')
addpath('models')
clc
volts2mu = 1;
TOL = 0.01;
trun = 800*40e-6;
ref_f = 2;
ref_0 = 0;
umax = 5;

matpath           = getMatPath();
dataroot          = fullfile(matpath, 'AFM_SS', 'System_Identification', 'data','data_xaxis'); 
expName           = '22-Jun-2016_exp01';
modFitName    = [expName, '.mat'];
modFitPath    = fullfile(dataroot, modFitName);
load(modFitPath, 'modelFit')


sys = ltiFit(modFitPath, 'SS02').sys;
Nd = 10;
sys_nodelay = sys;
sys_nodelay.InputDelay = 0;
sys.InputDelay = Nd;
sys = absorbDelay(sys);
Ts = sys.Ts;
PLANT = sys;

Ns = length(sys.b);

NsNd = Ns+Nd;



[uss_0, uss_f, ~, ~, xss]   = yss2uss(PLANT, ref_f, 0);
dcgain_sys = 1/(PLANT.c*xss);
x0 = xss*0;


% 3). LQR generation gain.        
% -------------------------------------------------------------------------
% -------------------- Constrained LQR Stuff ------------------------------
N_mpc = 8;
du_max   = 0.05;

% Pull out open-loop pole-zero information.
[wp_real_x, wz_real_x] = w_zp_real(sys);
rho_1 = wz_real_x(1)/wp_real_x(1);

% zeta_x = [.9, .8, .6, .5 .5];
zeta_x = [.9, .8, .7, .7 .7];
gams_x = [1.5, 1.5, 1.5, 1, 1];
rhos_x = [rho_1*1.0, 1, 1];
 
pint_x = 0.5;


%-----------------------------------------------------

if 1
    P_x    = getCharDes(sys_nodelay, gams_x, pint_x, zeta_x, rhos_x, .25);
    K_temp = place(sys_nodelay.a, sys_nodelay.b, P_x);
    [Q0, R1, K_lqr] = inverseLQR(sys_nodelay, K_temp);
    Q0 = blkdiag(Q0, zeros(Nd, Nd));
else
    P_x    = getCharDes(sys, gams_x, pint_x, zeta_x, rhos_x, .25);
    K_temp = place(sys.a, sys.b, P_x);
    [Q0, R1, K_lqr] = inverseLQR(sys, K_temp);
end

R1 = 1000;

% Q = blkdiag(Q_nodelay1, zeros(Nd,Nd));

sys_recyc = SSTools.deltaUkSys(sys);
Ns_mpc = size(sys_recyc.B, 1);
[Nx, Nu] = SSTools.getNxNu(sys_recyc);

Q1 = blkdiag(Q0, 0);
K_lqr = dlqr(sys_recyc.a, sys_recyc.b, Q1, R1);
Qp = dare(sys_recyc.a, sys_recyc.b, Q1, R1); 

R = 5;
mpc_on=0;



% *Optimal, open-loop* trajectory generation over all setpoints. This
% provides the basiline for comparison of the other methods. 

ref_s = [.4, .5];
N_traj = 600;
gamma = 100;
gam_s = linspace(gamma, 20000, 200);

% Form paramater classes
% 1. Time optimal

N_mpc_s = [4, 8, 12, 16, 20];
N_traj = 400;
trun = Ts*N_traj;
clear StepParamsMPC
clear StepParamsCLQR
clear build_clqr_trajs_obj
clear build_max_setpoints_obj
clear build_timeopt_trajs
clear StepData
clear StepParamsMPC

%%
clc
clear StepData
clear StepDataCLQR
clear StepDataTimeOpt
clear logger
clear ProgressBarFile
logfile = 'log-test.log';
logger = EchoFile(logfile);
ProgBar = @(max_iter, varargin)ProgressBarFile(max_iter, varargin{:}, 'logger', @logger.echo_file);
% Progbar(10)


step_params_timeopt = StepParamsTimeOpt(sys, ref_s, du_max, sys_nodelay, 10);
step_params_clqr  = StepParamsCLQR(sys_recyc, ref_s, du_max,Q1, gamma, PLANT, N_traj, 'condensed');

step_data_timeopt = StepDataTimeOpt(step_params_timeopt, 'savedata', true,...
    'file', 'data/timeopt_ref_data_test.mat', 'logger', @logger.echo_file,...
    'Progbar', ProgBar);
step_data_clqr = StepDataCLQR(step_params_clqr, 'savedata', true,...
    'file', 'data/clqr_ref_data_test.mat', 'logger', @logger.echo_file,...
    'Progbar', ProgBar);

%
% 1.----------- Generate for CLQR optimal trajectories ------------------


%%
clear StepDataCLQR
step_data_clqr = step_data_clqr.build_clqr_trajs('force', 1, 'verbose', 2);
%%
clear StepDataTimeOpt

clc
step_data_timeopt = step_data_timeopt.build_timeopt_trajs('force', 1, 'verbose', 2);

%%
clc
clear StepDataTimeOpt
clear StepData
figure
step_data_clqr.plot_single_traj(1)





%%
% expose variables for plotting.
clc
time_opt_settletime_s = step_data_timeopt.results.time_opt_settletime_s;
clqr_settletime_s = step_data_clqr.results.settle_times_opt_cell{1};



F300=figure(300); clf
colrs =  get(gca, 'colororder');
hands = [];
% nh = length(hands)
hands(1) = plot(ref_s, time_opt_settletime_s*1000, ':', 'LineWidth', 3);
set(hands(end), 'DisplayName', 'Time Optimal');
hold on
hands(2) = plot(ref_s, clqr_settletime_s*1000, 'Color', colrs(1, :), 'LineWidth', 2);
set(hands(end), 'DisplayName', 'CLQR')

legend(hands)
xlabel('reference [v]')
ylabel('settle-time [ms]')



