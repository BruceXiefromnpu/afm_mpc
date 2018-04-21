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



%%
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

trun = 400*Ts;

% mpcProb1 = sparseMPCprob(sys_recyc, N_mpc, Q1, Qp, R); 
mpcProb1 = condensedMPCprob(sys_recyc, N_mpc, Q1, Qp, R); 
mpcProb1.add_U_constraint('box', [-du_max, du_max]);



sim_struct = struct('PLANT', PLANT, 'trun', trun, 'mpcProb1', mpcProb1,...
                    'du_max', du_max, 'mpc_on', mpc_on,...
                    'xss', xss, 'Nx', Nx);

ref_range = [0.01, 15];
ref_step = .1;
ref_s = ref_range(1):ref_step:ref_range(2);

% *Optimal, open-loop* trajectory generation over all setpoints. This
% provides the basiline for comparison of the other methods. 
N_traj = 400;
gamma = 100;
gam_s = linspace(100, 20000, 5);
max_setpoints_linear = 0*gam_s;



%%
data_struct.N_mpc_s = [4, 8];
data_struct.sys = sys_recyc;
data_struct.Q   = Q1;
data_struct.ref_s = ref_s;
data_struct.gam_s = gam_s(1:3);
data_struct.sim_struct = sim_struct;
data_struct.mpc_on    = 0;
data_struct.verbose   = 1;
data_struct.savedata  = 1;
data_struct.file      = 'data/max_ref_data_dalay.mat';
data_struct.fig_files = ['figures/max_ref_dataf1.fig', 'figures/max_ref_data2.fig'];


clear build_clqr_trajs
clear opt_traj_jen
data_struct_clqr = data_struct;
data_struct_clqr.gam_s_clqr = 100;
data_struct_clqr.file = 'data/clqr_data_tmp.mat';
data_struct_clqr.du_max = du_max;
data_struct_clqr.N_traj = 600;
data_struct_clqr.mpc_mode = 'condensed';

clqr_data = build_clqr_trajs(data_struct_clqr, 'force', 0);

F200 = figure(200); hold on;
settle_times_clqr = clqr_data.settle_times_opt_save{1};
plot(ref_s, settle_times_clqr)

% settle_times_opt_save = cell(1, length(gam_s));
% opt_trajs_save = cell(1, length(gam_s));
% hands = [];
% F200 = figure(200); hold on
% colrs = get(gca, 'colororder');
% load_opt = 0;
% for iter = 1:length(gam_s)
%     gamma = gam_s(iter);
%     if load_opt
%         load('clqr_opt_data.mat')
%         opt_trajs_save = clqr_opt_data.opt_trajs_save;
%         settle_times_opt_save = clqr_opt_data.settle_times_opt_save;
%     else
%         [traj_s, settle_times_opt] = opt_traj_gen(Q1, gamma, N_traj, sys_recyc, ...
%                                 ref_s, 0, 'uMax', du_max, 'verbose', 0);
%         settle_times_opt_save{iter} = settle_times_opt;
%         opt_trajs_save{iter} = traj_s;
%     end
%     change_current_figure(F200); 
%     hands(iter) = plot(ref_s, settle_times_opt_save{iter}*1000, 'LineWidth', 2);
% 
%     set(hands(iter), 'DisplayName', sprintf('$\\gamma = %.0f$', gamma),...
%         'Color', colrs(iter, :));
%     leg = legend(hands);
%     set(leg, 'interpreter', 'latex');
%     ylabel('settle time [ms]', 'FontSize', 16)
%     xlabel('setpoint', 'FontSize', 16)
%     drawnow()
%     grid on
% end
% leg.FontSize = 16;
% leg.Location = 'Northwest';

%% for all the ref_f_s, get time-optimal settle times
g_eject = zpk(exp(-wp_real_x(1)*Ts), exp(-wz_real_x(1)*Ts), 1, sys_nodelay.Ts);
g_eject = g_eject/dcgain(g_eject);
sys_eject = minreal(sys_nodelay*g_eject);
sys_recyc_eject = SSTools.deltaUkSys(sys_eject);
Nx_eject = SSTools.getNxNu(sys_recyc_eject);
x0_eject = Nx_eject*0;
time_opt_settletime_s = ref_s*0;
warning('off', 'MATLAB:nargchk:deprecated')

clear timeOptBisect
toBisect = TimeOptBisect(sys_recyc_eject, du_max);
toBisect.maxIter = 50;
vb = 1;
if vb > 1
    F1000 = figure(1000)
end
for k=52:length(ref_s)
    xf = Nx_eject*ref_s(k);
    [xx, uu, S]=toBisect.time_opt_bisect(x0_eject, xf);
    time_opt_settletime_s(k) = uu.Time(end) + 10*sys_recyc_eject.Ts;
    
    if vb > 0
 
        if vb > 1 
           if mod(k, 20)==0
            for i=1:k
                close(figure(1000+i))
            end
        end
            figure(1000+k)
        else
            change_current_figure(F1000);
        end
        set(gcf, 'Position', [-981 49 894 976])
        % Extend the trajectory, so we can see if its actually at steady state;
        t_ext = [uu.Time(end)+Ts:Ts:uu.Time(end)*2]';
        t = [uu.Time; t_ext];
        u = [uu.Data; 0*t_ext];
        y = lsim(sys_recyc_eject, u, t, zeros(12,1));
        fprintf('optimal N=%.0f\n', length(uu.Data))
        subplot(311); 
        plot(t, y)
        hold on
        xlm = xlim;
        plot(xlm, [ref_s(k), ref_s(k)], '--')

        subplot(312)
        plot(t, u)
        hold on;
        
        subplot(313)
        plot(t, cumsum(u))
        hold on
        
        drawnow()
    end
%     keyboard
end
%%
figure(200)

% nh = length(hands)
hands(end+1) = plot(ref_s, time_opt_settletime_s*1000, ':', 'LineWidth', 3)
set(hands(end), 'DisplayName', 'Time Optimal')
legend(hands)

%%
data_struct.N_mpc_s = [4, 8];
data_struct.sys = sys_recyc;
data_struct.Q   = Q1;
data_struct.ref_s = ref_s;
data_struct.gam_s = gam_s(1:3);
data_struct.sim_struct = sim_struct;
data_struct.mpc_on = 0;
data_struct.verbose = 1;
data_struct.savedata = 1;


data_struct.file = 'data/max_ref_data_dalay.mat';
data_struct.file = 'data/max_ref_data_dalay.mat';
data_struct.fig_files = ['figures/max_ref_dataf1.fig', 'figures/max_ref_data2.fig'];
clc
clear build_max_setpoints
build_max_setpoints(data_struct)
%%
saveon = 1;
clqr_opt_data.opt_trajs_save = opt_trajs_save;
clqr_opt_data.settle_times_opt_save = settle_times_opt_save;
if saveon
    saveas(F200, 'figures/opttraj_setpoint_vs_ts.svg')
%     save('data/clqr_opt_data.mat', 'clqr_opt_data');
end


%%
clc
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

