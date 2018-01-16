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

fname_lin = 'data/max_sp_data_lin_CCTA.mat';
fname_mpc = 'data/max_sp_data_mpc_CCTA.mat';
fname_clqr = 'data/clqr_ref_data_CCTA.mat';

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
du_max   = 0.1;

% Pull out open-loop pole-zero information.
[wp_real_x, wz_real_x] = w_zp_real(sys);
rho_1 = wz_real_x(1)/wp_real_x(1);


% zeta_x = [.9, .8, .7, .7 .7]; original 
% gams_x = [1.5, 1.5, 1.5, 1, 1]; original
zeta_x = [.9, .8, .7, .7 .5];
gams_x = [1.5, 1.5, 1.5, 1.25, 1];
rhos_x = [rho_1*1.0, 1, 1];
 
pint_x = 0.5;


%-----------------------------------------------------


P_x    = getCharDes(sys_nodelay, gams_x, pint_x, zeta_x, rhos_x, .25);
K_temp = place(sys_nodelay.a, sys_nodelay.b, P_x);
[Q0, R1, K_lqr] = inverseLQR(sys_nodelay, K_temp);
Q0 = blkdiag(Q0, zeros(Nd, Nd));



R1 = 1000;

% Q = blkdiag(Q_nodelay1, zeros(Nd,Nd));

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

clc

R = 5;
mpc_on=0;

ref_range = [0.01, 15];
ref_step = .1;
ref_s = ref_range(1):ref_step:ref_range(2);

% *Optimal, open-loop* trajectory generation over all setpoints. This
% provides the basiline for comparison of the other methods. 
% ref_s = [1.5];
N_traj = 600;
gamma = 100;
gam_s = linspace(gamma, 20000, 200); % original
% gam_s = [1, 100, 1000, 2500, 5000, 10000];

% Form paramater classes
% 1. Time optimal

N_mpc_s = [4, 8, 12, 16, 20]; % original 
% N_mpc_s = [12, 18, 24];
N_traj = 800;
trun = Ts*N_traj;

step_params_lin = StepParamsLin(sys_recyc, ref_s, du_max,Q1, gam_s, PLANT, trun);
step_params_mpc = StepParamsMPC(sys_recyc, ref_s, du_max,Q1, gam_s, PLANT, N_mpc_s, 'condensed', trun);
max_sp_data_lin = StepData(step_params_lin, 'savedata', true,...
    'file', fname_lin);
max_sp_data_mpc = StepData(step_params_mpc, 'savedata', true,...
    'file', fname_mpc);

% step_params_timeopt = StepParamsTimeOpt(sys, ref_s, du_max, sys_nodelay, 10);
step_params_clqr  = StepParamsCLQR(sys_recyc, ref_s, du_max,Q1, gamma, PLANT, N_traj, 'condensed');
step_data_clqr = StepDataCLQR(step_params_clqr, 'savedata', true,...
    'file', fname_clqr);


%
% 1.----------- Generate for CLQR optimal trajectories ------------------

logfile = 'log_clqr.log';
echo_file([], logfile);
logstart=2;


%%
% 1.----------- Generate LIN max setpoints --------------------------------
tic
try
    max_sp_data_lin = build_max_setpoints(max_sp_data_lin, 'fid', fid,'verbose', 1, 'force', 0);
    fprintf(fid, 'Finished building max setpoints, linear. Total time = %.2f\n\n', toc);
catch ME
    errMsg = getReport(ME, 'extended', 'hyperlinks', 'off');
     fprintf(fid, 'Failed to build max setpoints, linear: \n\n%s', errMsg);
end

% 2.----------- Generate  mpc max setpoints -------------------------------

tic
try
    max_sp_data_mpc = build_max_setpoints(max_sp_data_mpc, 'fid', fid, 'verbose', 1, 'force', 1);
    fprintf(fid, 'Finished building max setpoints, mpc. Total time = %.2f\n\n', toc);
catch ME
    errMsg = getReport(ME,  'extended','hyperlinks', 'off');
    fprintf(fid, 'Failed to build max setpoints, mpc: %s\n\n', errMsg);
end

% 4.----------- Generate CLQR  trajectories -------------------------------
tic
try
    step_data_clqr = build_clqr_trajs(step_data_clqr, 'force', 1, 'fid', fid, 'verbose', 1);
    fprintf(fid, 'Finished building clqr_data. Total time = %.2f\n', toc);
catch ME
    errMsg = getReport(ME, 'extended', 'hyperlinks', 'off');
    fprintf(fid, 'Failed to build clqr_data: \n%s', errMsg);
end




% if fid >1
% %     fclose(fid);
%     to = { 'abraker@fastmail.com', 'robr9299@colorado.edu'};
%     ssendmail('MATLAB Report:compare_maximum_ref.mm', logfile, 'to', to)
% %     ssendmail('test from matlab', 'body of test from matlab', 'to',...
% %         'robr9299@colorado.edu', 'attachments', {'figures/cp_traj.svg'})
% end



% expose variables for plotting.
try
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

    saveon = 1;

    if saveon
        saveas(F300, 'figures/clqrTimeOpt_sp_vs_ts_CCTA.svg')
    end
    
    % Plot maximum reference vs gamma
    F10 = figure(10); clf; hold on;
    colrs =  get(gca, 'colororder');
    hands_mxsp = []
    hands_mxsp(1) = plot(max_sp_data_lin.params.gam_s, max_sp_data_lin.results{1}.max_setpoints,...
                         '-', 'Color', colrs(1, :),...
                'LineWidth', 2);
    set(hands_mxsp(1), 'DisplayName', 'LQR Lin + Sat');        
    for k = 1:length(max_sp_data_mpc.params.N_mpc_s)

        hands_mxsp(k+1) = plot(max_sp_data_mpc.params.gam_s, max_sp_data_mpc.results{k}.max_setpoints,...
                           '-', 'Color', colrs(k+1, :), 'LineWidth', 2);

        stit = sprintf('MPC, N=%.0f', max_sp_data_mpc.params.N_mpc_s(k));
        set(hands_mxsp(k+1), 'DisplayName', stit);
    end
    legend(hands_mxsp)

    xlabel('$\gamma$', 'interpreter', 'latex', 'FontSize', 14);
    ylabel('Max Ref', 'interpreter', 'latex', 'FontSize', 14);

    saveas(F10, 'figures/max_ref_vs_gamma_lin_MPC_CCTA.svg')

catch
    fprintf(fid, 'Error making plots\n')
end

% ======================================================================= %
%                                                                         %
% ------------------Settle-By-Maximum-Acheivable-Reference -------------- %
%                                                                         %
% ======================================================================= %
clc
clear TsByMaxRefParams
logfile2 = 'log-tsbyrefMax.txt';
fid = fopen('log-tsbyrefMax.txt', 'w');

% fid = 1;
rmax_s = [1, 1.5, 2.5, 5.0, 10];
ts_by_rmax_lin = TsByMaxRefParams(max_sp_data_lin, rmax_s, 'data/ts_byrefmax_lin.mat');
ts_by_rmax_mpc = TsByMaxRefParams(max_sp_data_mpc, rmax_s, 'data/ts_byrefmax_mpc.mat');
tic
try
    ts_by_rmax_lin = ts_by_rmax_lin.run_ts_by_refs('verbose', 1, 'fid', fid);
    fprintf(fid, 'Finished running Ts by Max Ref, linear. Total time = %.2f\n\n', toc);
catch ME
    errMsg = getReport(ME,  'extended','hyperlinks', 'off');
    fprintf(fid, 'Failed to run Ts by max setpoints, linear: %s\n\n', errMsg);
end

tic
try
    ts_by_rmax_mpc = ts_by_rmax_mpc.run_ts_by_refs_mpc('verbose', 1, 'fid', fid);
    fprintf(fid, 'Finished running Ts by Max Ref, MPC. Total time = %.2f\n\n', toc);
catch ME
    errMsg = getReport(ME,  'extended','hyperlinks', 'off');
    fprintf(fid, 'Failed to run Ts by max setpoints, mpc: %s\n\n', errMsg);
end



% if fid >1
% %     fclose(fid);
%     to = { 'abraker@fastmail.com', 'robr9299@colorado.edu'};
%     ssendmail('MATLAB Report:compare_maximum_ref.mm', logfile2, 'to', to)
% %     ssendmail('test from matlab', 'body of test from matlab', 'to',...
% %         'robr9299@colorado.edu', 'attachments', {'figures/cp_traj.svg'})
% end
% 


try
    clc
    clear TsByMaxRefParams
    figure(20)

    [ax, hands, leg] = ts_by_rmax_lin.plot_ts_v_r2max(gca);
    

    % Copy the ref-vs-ts figure;
    f400 = copyobj(F300, 0);
    hold on

    % hands2 = f400.Children.Children;
    hands2 = f400.CurrentAxes.Children;
    hands3 = gobjects(length(ts_by_rmax_lin.rmax_s), 1)

    for k=1:length(ts_by_rmax_lin.rmax_s)
        rmax = ts_by_rmax_lin.rmax_s(k);
        settle_times = ts_by_rmax_lin.results{k}.t_settle_s;
        ref_s = ts_by_rmax_lin.results{k}.ref_s_to_rmax;
        gamma = ts_by_rmax_lin.results{k}.gamma;

        stit = sprintf('LIN, $r_{max}=%.2f$, $\\gamma=%.0f$', rmax, gamma);
        hands3(k) = plot(ref_s, settle_times*1000);
        set(hands3(k), 'DisplayName', stit, 'LineWidth', 2);
    end
    
    hleg = legend([hands2; hands3]);
    set(hleg, 'interpreter', 'latex')


    clc
    num_N_mpcs = length(ts_by_rmax_mpc.StepData.params.N_mpc_s);
    fighands = gobjects(num_N_mpcs, 1);

    for mpc_iter=1:num_N_mpcs
        fighands(mpc_iter) = copyobj(F300, 0);
        hold on

        % hands2 = f400.Children.Children;
        hands2 = fighands(mpc_iter).CurrentAxes.Children;
        hands3 = gobjects(length(ts_by_rmax_mpc.rmax_s), 1)


        N_mpc = ts_by_rmax_mpc.StepData.params.N_mpc_s(mpc_iter)

        for k=1:length(ts_by_rmax_mpc.rmax_s)
            rmax = ts_by_rmax_mpc.rmax_s(k);
            settle_times = ts_by_rmax_mpc.results{mpc_iter}{k}.t_settle_s;
            ref_s = ts_by_rmax_mpc.results{mpc_iter}{k}.ref_s_to_rmax;
            gamma = ts_by_rmax_mpc.results{mpc_iter}{k}.gamma;

            stit = sprintf('MPC: N=%.0f, $r_{max}=%.2f$, $\\gamma=%.0f$', N_mpc, rmax, gamma);
            hands3(k) = plot(ref_s, settle_times*1000);
            set(hands3(k), 'DisplayName', stit, 'LineWidth', 2);
        end

        title(sprintf('Nmpc = %.0f', N_mpc))


        hleg = legend([hands2; hands3]);
        set(hleg, 'interpreter', 'latex')

    end

catch
    fprintf(fid, 'error making Max ref figures\n')
end

fclose(fid)



