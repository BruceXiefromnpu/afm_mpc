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

% mpcProb0 = condensedMPCprob(sys_recyc, N_mpc, Q1, Qp, R1);
% mpcProb0.Ainq = [eye(N_mpc); -eye(N_mpc)];
% mpcProb0.binq = [zeros(2*N_mpc, 1)+du_max];




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

% ref_range = [0.01, 15];
% ref_step = .1;
% ref_s = ref_range(1):ref_step:ref_range(2);

% *Optimal, open-loop* trajectory generation over all setpoints. This
% provides the basiline for comparison of the other methods. 
ref_s = [0.1:0.01:15];
N_traj = 600;
gamma = 100;
gam_s = linspace(gamma, 20000, 200);

% Form paramater classes
% 1. Time optimal

N_mpc_s = [4, 8, 12, 16, 20];
N_traj = 800;
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

step_params_timeopt = StepParamsTimeOpt(sys, ref_s, du_max, sys_nodelay, 10);
step_params_lin   = StepParamsLin(sys_recyc, ref_s, du_max,Q1, gam_s, PLANT, trun);
step_params_mpc   = StepParamsMPC(sys_recyc, ref_s, du_max,Q1, gam_s, PLANT, N_mpc_s, 'condensed', trun);
step_params_clqr  = StepParamsCLQR(sys_recyc, ref_s, du_max,Q1, gamma, PLANT, N_traj, 'condensed');

step_data_timeopt = StepData(step_params_timeopt, 'savedata', true,...
    'file', 'data/timeopt_ref_data.mat');
step_data_clqr = StepDataCLQR(step_params_clqr, 'savedata', true,...
    'file', 'data/clqr_ref_data.mat');
max_sp_data_lin = StepData(step_params_lin, 'savedata', false,...
    'file', 'data/max_sp_data_lin.mat');
max_sp_data_mpc = StepData(step_params_mpc, 'savedata', false,...
    'file', 'data/max_sp_data_mpc.mat');

%
% 1.----------- Generate for CLQR optimal trajectories ------------------

% fid = fopen('log2.txt', 'w+');
fid = 1;
%%
tic
try
    step_data_clqr = build_clqr_trajs(step_data_clqr, 'force', 0, 'fid', fid, 'verbose', 2);
    fprintf(fid, 'Finished building clqr_data. Total time = %.2f\n', toc);
catch ME
    errMsg = getReport(ME, 'extended', 'hyperlinks', 'off');
    fprintf(fid, 'Failed to build clqr_data: \n%s', errMsg);
end
%%
% 2.----------- Generate LIN max setpoints --------------------------
tic
try
%     clear find_ref_max
%     clear build_max_setpoints
%     max_sp_data_lin.params.ref_s = linspace(1, 5, 10);
%     max_sp_data_lin.params.gam_s = linspace(100, 10000, 10);
    max_sp_data_lin = build_max_setpoints(max_sp_data_lin, 'fid', fid, 'force', 0);
    fprintf(fid, 'Finished building max setpoints, linear. Total time = %.2f\n\n', toc);
catch ME
    errMsg = getReport(ME, 'extended', 'hyperlinks', 'off');
     fprintf(fid, 'Failed to build max setpoints, linear: \n\n%s', errMsg);
end
%%
% 3.----------- Generate  mpc max setpoints --------------------------

tic
try
    max_sp_data_mpc = build_max_setpoints(max_sp_data_mpc, 'fid', fid, 'verbose', 2);
    fprintf(fid, 'Finished building max setpoints, mpc. Total time = %.2f\n\n', toc);
catch ME
    errMsg = getReport(ME,  'extended','hyperlinks', 'off');
    fprintf(fid, 'Failed to build max setpoints, mpc: %s\n\n', errMsg);
end
%%
% 4.----------- Generate for time-optimal trajectories ------------------
clear build_timeopt_trajs
try
%    step_data_timeopt.params.ref_s = ref_s(50:52)
   [step_data_timeopt, status] = build_timeopt_trajs(step_data_timeopt, 'force', 0, 'fid', fid, 'max_iter', 50);
   if status
       fprintf(fid, 'build_timeopt_trajs failed with exit status %0.0f\n', status);
   else
       fprintf(fid, 'Finished building time-optimal trajectories. Total time = %.2f\n', toc);
   end
   catch ME
    exception = MException.last;
    errMsg = getReport(ME, 'extended', 'hyperlinks', 'off');
    fprintf(fid, 'Failed to build time-optimal trajectories:\n%s', errMsg);
end

if fid >1
    fclose(fid);
    msg = freadf('log.txt');
    to = { 'abraker@fastmail.com', 'robr9299@colorado.edu'};
    ssendmail('MATLAB Report:compare_maximum_ref.mm', msg, 'to', to)
%     ssendmail('test from matlab', 'body of test from matlab', 'to',...
%         'robr9299@colorado.edu', 'attachments', {'figures/cp_traj.svg'})
end


%%
% expose variables for plotting.
% time_opt_settletime_s = step_data_timeopt.results.time_opt_settletime_s;
clqr_settletime_s = step_data_clqr.results.settle_times_opt_cell{1};



figure(300); clf
colrs =  get(gca, 'colororder');
hands = []
% nh = length(hands)
% hands(1) = plot(ref_s, time_opt_settletime_s*1000, ':', 'LineWidth', 3);
% set(hands(end), 'DisplayName', 'Time Optimal')
hold on
hands(2) = plot(ref_s, clqr_settletime_s*1000, 'Color', colrs(1, :), 'LineWidth', 2);
set(hands(end), 'DisplayName', 'CLQR')

legend(hands)

%%
saveon = 1;
clqr_opt_data.opt_trajs_save = opt_trajs_save;
clqr_opt_data.settle_times_opt_save = settle_times_opt_save;
if saveon
    saveas(F200, 'figures/opttraj_setpoint_vs_ts.svg')
%     save('data/clqr_opt_data.mat', 'clqr_opt_data');
end
%%
% Plot maximum reference vs gamma
F10 = figure(10); clf; hold on;
colrs =  get(gca, 'colororder');
hands_mxsp = []
hands_mxsp(1) = plot(max_sp_data_lin.params.gam_s, max_sp_data_lin.results.max_setpoints,...
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

saveas(F10, 'figures/max_ref_vs_gamma_lin_MPC.svg')
%%
clc
clear TsByMaxRefParams
fid = fopen('log-tsbyrefMax.txt', 'w+');
% fid = 1;
rmax_s = [1, 2.1, 4.8, 5.5, 10];
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
fprintf(fid, 'Finished running Ts by Max Ref, MPC. Total time = %.2f\n\n', toc);
ts_by_rmax_mpc = ts_by_rmax_lin.run_ts_by_refs_mpc('verbose', 1, 'fid', fid);

catch ME
    errMsg = getReport(ME,  'extended','hyperlinks', 'off');
    fprintf(fid, 'Failed to run Ts by max setpoints, mpc: %s\n\n', errMsg);
    
end

if fid >1
    fclose(fid);
    msg = freadf('log-tsbyrefMax.txt');
    to = { 'abraker@fastmail.com', 'robr9299@colorado.edu'};
    ssendmail('MATLAB Report:compare_maximum_ref.mm', msg, 'to', to)
%     ssendmail('test from matlab', 'body of test from matlab', 'to',...
%         'robr9299@colorado.edu', 'attachments', {'figures/cp_traj.svg'})
end



%%
clc
clear TsByMaxRefParams
figure(20)

[ax, hands, leg] = ts_by_rmax_lin.plot_ts_v_r2max(gca);
%%
% Now load in the data generated by build_max_setpoints.m
load('data/max_ref_data_dalay.mat')
gam_s = max_ref_data.gam_s;
max_setpoints_linear = max_ref_data.max_setpoints_lin;
% Load CLQR opt data.
load('clqr_opt_data.mat')
ts_opt = clqr_opt_data.settle_times_opt_save{1};


% Now, lets ask the question: Take a set of LINEAR LQR controllers
% parameterized by gamma. Then for each K(gamma), we now have a maximum
% setpoint value we can visit. So for each one of these K(gamma)'s, lets
% plot the settling time vs setpoint size. I expect to find that settling
% time decreases as max-setpoint decreases. 

F400 = figure(400); clf
ylabel('Settle-time perc. increase', 'FontSize', 16)
xlabel('setpoint', 'FontSize', 16)
hands_perc = [];
grid on

F300 = figure(300); clf
hands = [];
hands(1) = plot(ref_s, settle_times_opt_save{1}*1000, 'LineWidth', 2,...
                'Color', colrs(1, :));
set(hands(1), 'DisplayName', 'CLQR opt')

ylabel('settle time [ms]', 'FontSize', 16)
xlabel('setpoint', 'FontSize', 16)
grid







rmax_s = [1, 2.1, 4.8];
subs = [221, 222, 223, 224];
% figure(100)
for rmax_iter = 1:length(rmax_s)
% 1. linear guy
%     rmax = 1.0;
    rmax = rmax_s(rmax_iter);
    
    kk = find(max_setpoints_linear >=rmax, 1, 'first');

    gamma = gam_s(kk);
    % 
    % sim_struct = struct('PLANT', PLANT, 'trun', trun, 'mpcProb1', mpcProb1,...
    %                     'du_max', du_max, 'mpc_on', mpc_on,...
    %                     'xss', xss, 'Nx', Nx);
    K_lqr = dlqr(sys_recyc.a, sys_recyc.b, Q1, gamma);

    t_settle_lin_s = zeros(1, kk);
    k_ref = find(ref_s == max_setpoints_linear(kk), 1, 'first');
    labs = [];

    for k=1:k_ref
        mpc_on = 0;
        ref_f = ref_s(k);
        sim('MPC_fp')

        y1 = y_mpcDist; 
        [t_settle, k_s] = settle_time(y1.time, y1.Data, ref_f, 0.01*ref_f,...
                                      [], [], 30);
        t_settle_lin_s(k) = t_settle;
        
        if k == 22 && rmax_iter == 3
            figure(500); 
            set(gcf, 'Position', [-971 66 876 939]);
            
            subplot(321);
            plot(y1.time, y1.Data, 'Color', colrs(rmax_iter+1, :),...
                'LineWidth', 2)
            hold on
            opt = opt_trajs_save{1};
            yopt = opt.Y_vec_s;
            plot(yopt{k}.Time, yopt{k}.Data, '--', 'Color', colrs(1,:),...
                'LineWidth', 2)
            ylabel('y(t)')
            % U
            subplot(323)
            plot(u_mpcDist.Time, u_mpcDist.Data, 'Color', colrs(rmax_iter+1, :),...
                'LineWidth', 2)
            hold on
            uopt = opt.U_vec_s{k};
            uopt.Data = cumsum(uopt.Data);
            plot(uopt.Time, uopt.Data, '--', 'Color', colrs(1,:),...
                'LineWidth', 2)
            xlabel('time')
            ylabel('u')
            
            % DU
            subplot(325)
            plot(du_mpcDist.Time, du_mpcDist.Data, 'Color', colrs(rmax_iter+1, :),...
                'LineWidth', 2)
            hold on
            duopt = opt.U_vec_s{k};
            plot(duopt.Time, duopt.Data, '--', 'Color', colrs(1,:),...
                'LineWidth', 2)
            xlabel('time')
            ylabel('du')

        elseif  k == 30 && rmax_iter == 3
            
            figure(500); 
            subplot(322)
            plot(y1.time, y1.Data, 'Color', colrs(rmax_iter+1, :),...
                'LineWidth', 2)
            hold on
             opt = opt_trajs_save{1};
            yopt = opt.Y_vec_s;
            plot(yopt{k}.Time, yopt{k}.Data, '--', 'Color', colrs(1,:),...
                'LineWidth', 2)
            ylabel('y(t)')
            
            % U
            subplot(324)
            plot(u_mpcDist.Time, u_mpcDist.Data, 'Color', colrs(rmax_iter+1, :),...
                'LineWidth', 2)
            hold on
            uopt = opt.U_vec_s{k};
            uopt.Data = cumsum(uopt.Data);
            plot(uopt.Time, uopt.Data, '--', 'Color', colrs(1,:),...
                'LineWidth', 2)
            xlabel('time')
            ylabel('u(t)')

            % DU
            subplot(326)
            plot(du_mpcDist.Time, du_mpcDist.Data, 'Color', colrs(rmax_iter+1, :),...
                'LineWidth', 2)
            hold on
            duopt = opt.U_vec_s{k};
            plot(duopt.Time, duopt.Data, '--', 'Color', colrs(1,:),...
                'LineWidth', 2)
            xlabel('time')
            ylabel('du')
            
        end
            
%         figure(101); hold on
%         plot(y1.Time, y1.Data)
    end

    clc
    figure(F300)
%     subplot(subs(rmax_iter)); 
    hold on
    hands(rmax_iter+1) = plot(ref_s(1:k_ref), t_settle_lin_s*1000,...
                    'LineWidth', 2, 'Color', colrs(rmax_iter+1, :));

    leg_str = sprintf('linear: $\\gamma=%.0f$, $r_{max}=%.1f$', gamma, rmax);
    set(hands(rmax_iter+1), 'DisplayName', leg_str);
    leg = legend(hands);
    set(leg, 'interpreter', 'latex');
    leg.FontSize = 14;
    leg.Position = [0.4842 0.1500 0.4845 0.2281];
    drawnow()
    
    
    % Plot percentage increase/decrease.
        
    ts_opt_k = ts_opt(1:length(t_settle_lin_s));
    perc_inc = t_settle_lin_s./ts_opt_k;
    
    change_current_figure(F400);
    hold on;
    hands_perc(rmax_iter) = plot(ref_s(1:k_ref), perc_inc,...
                        'LineWidth', 2, 'Color', colrs(rmax_iter+1, :));
    leg_str_perc = sprintf('$\\gamma = %.0f$, $r_{max} = %.1f$', gamma, rmax);
    set(hands_perc(rmax_iter), 'DisplayName', leg_str_perc);
    
    leg_perc = legend(hands_perc);
    set(leg_perc, 'interpreter', 'latex', 'FontSize', 14);
    
%     grid on
%     ylabel('t-settle')
%     xlabel('ref')

    % Now do the MPC
%     t_settle_s = zeros(length(N_mpc_s), kk);
%     mpc_on = 1;
%     for mpc_iter = 1:length(N_mpc_s(1:end-1))
%         N_mpc = N_mpc_s(mpc_iter);
%         kk = find(max_setpoints_mpc(mpc_iter, :) >=rmax, 1, 'first');
%         gamma = gam_s(kk+1);
% 
%         Qp = dare(sys_recyc.a, sys_recyc.b, Q1, gamma); 
%         % sim_struct.mpcProb1 = sparseMPCprob(sys_recyc, N_mpc, Q1, Qp, gamma);
%         mpcProb1 = condensedMPCprob(sys_recyc, N_mpc, Q1, Qp, gamma);
%         mpcProb1.add_U_constraint('box', [-du_max, du_max]);
% 
%         K_lqr = dlqr(sys_recyc.a, sys_recyc.b, Q1, gamma)*0;
% 
% 
%         t_settle_s = zeros(1, kk);
%         k_ref = find(ref_s >=rmax, 1, 'first');
%         for k=1:k_ref
% 
%             ref_f = ref_s(k)
%             sim('MPC_fp')
% 
%             y1 = y_mpcDist; 
%             [t_settle, k_s] = settle_time(y1.time, y1.Data, ref_f, 0.01*ref_f,...
%                                           [], [], 30);
% 
%             t_settle_s(mpc_iter, k) = t_settle;
%             figure(mpc_iter); hold on;
%             plot(y1.Time, y1.Data)
%             drawnow()
%         end
%         
%         figure(100)
%         subplot(subs(rmax_iter));
%         hold on
%         
%         title(sprintf('rmax = %.1f', rmax));
%         labs(mpc_iter+1) = plot(ref_s(1:k_ref), t_settle_s(mpc_iter,:));
%         legstr = sprintf('Mpc: N= %.0f, gam=%.0f', N_mpc, gamma);
%         set(labs(mpc_iter+1), 'DisplayName', legstr)
%         legend(labs)
%     end

end
%%

if saveon

   saveas(gcf, 'figures/cp_traj.svg')
end
%%
if saveon
    saveas(F300, 'figures/cp_clqropt_vs_linmaxset_TS-s.svg')
    saveas(F400, 'figures/perc_increase_lin_over_clqr.svg');
end

