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
expName           = ['22-Jun-2016_exp01'];
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

mpcProb0 = condensedMPCprob(sys_recyc, N_mpc, Q1, Qp, R1);
mpcProb0.Ainq = [eye(N_mpc); -eye(N_mpc)];
mpcProb0.binq = [zeros(2*N_mpc, 1)+du_max];

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

gam_s = linspace(100, 20000, 50);
gam_s = gam_s(1:25);

R = 5;

mpc_on=0;

trun = 400*Ts;

% mpcProb1 = sparseMPCprob(sys_recyc, N_mpc, Q1, Qp, R); 
mpcProb1 = condensedMPCprob(sys_recyc, N_mpc, Q1, Qp, R); 
mpcProb1.add_U_constraint('box', [-du_max, du_max]);

max_setpoints_linear = 0*gam_s;

sim_struct = struct('PLANT', PLANT, 'trun', trun, 'mpcProb1', mpcProb1,...
                    'du_max', du_max, 'mpc_on', mpc_on,...
                    'xss', xss, 'Nx', Nx);

ref_range = [0.01, 5];
ref_step = .1;
ref_s = ref_range(1):ref_step:ref_range(2);




upd = progressbar(length(gam_s));
F1 = figure(1); clf;
F2 = figure(2); clf; hold on;
colrs = get(gca, 'colororder');
ylabel('settle time [ms]')
xlabel('ref_f')
reverseStr = '';


for gam_iter = 7:7
%     1:length(gam_s)
    gamma = gam_s(gam_iter);
    sim_struct.K_lqr = dlqr(sys_recyc.a, sys_recyc.b, Q1, gamma);
    [ref_max, t_settle_s] = find_ref_max(sim_struct, ref_s,...
                            'verbose', 1, 'fig_base', 10*gam_iter);
    
    change_current_figure(F2);
    k_max = find(t_settle_s == 0, 1, 'first') -1;
    plot(ref_s(1:k_max), t_settle_s(1:k_max)*1e3)
    ylim([0, 10])
    max_setpoints_linear(gam_iter) = ref_max;

    if exist('hlin', 'var')
        delete(hlin)
        drawnow()
    end

    change_current_figure(F1)
    hlin = plot(gam_s(1:gam_iter), max_setpoints_linear(1:gam_iter),...
        '-o', 'Color', colrs(1,:), 'LineWidth', 2);
    drawnow()
    hold on

    upd(gam_iter);
%     keyboard
end

hlin.DisplayName = 'LQR lin + sat';

%%
% ----------------------- Same thing, but for MPC ----------------------- %

sim_struct.mpc_on = 1;
% upd = textprogressbar(length(gam_s));

% N_mpc = 8;
N_mpc_s = [4, 8, 12, 16, 20];
max_setpoints_mpc = repmat(0*gam_s, length(N_mpc_s), 1);


F20 = figure(20); clf; hold on;
colrs = get(gca, 'colororder');
ylabel('settle time [ms]')
xlabel('ref_f')


hmpc_s = [];
for mpc_iter = 1:length(N_mpc_s)
    N_mpc = N_mpc_s(mpc_iter);
    colr = colrs(mpc_iter+1, :);
    start_str = sprintf('N = %.0f', N_mpc);
    upd = progressbar(length(gam_s), 'bar_len', 25, 'start_str', start_str);
    
    % Inititalize the progress bar.
    upd(0)
    for gam_iter = 1:length(gam_s)
        gamma = gam_s(gam_iter);

        % -------------- Update the MPC parameters ---------------------- %
        % must re-compute terminal cost when gamma changes!
        Qp = dare(sys_recyc.a, sys_recyc.b, Q1, gamma); 
        % sim_struct.mpcProb1 = sparseMPCprob(sys_recyc, N_mpc, Q1, Qp, gamma);
        sim_struct.mpcProb1 = condensedMPCprob(sys_recyc, N_mpc, Q1, Qp, gamma);
        sim_struct.mpcProb1.add_U_constraint('box', [-du_max, du_max]);

        sim_struct.K_lqr = dlqr(sys_recyc.a, sys_recyc.b, Q1, gamma)*0;

        [ref_max, t_settle_s] = find_ref_max(sim_struct, ref_s,...
                            'verbose', 1, 'fig_base', 10*gam_iter);

        
        change_current_figure(F20);
        k_max = find(t_settle_s == 0, 1, 'first') -1;
        plot(ref_s(1:k_max), t_settle_s(1:k_max)*1e3)
        ylim([0, 10])
        
        max_setpoints_mpc(mpc_iter, gam_iter) = ref_max;
        
        if length(hmpc_s) == mpc_iter
            delete(hmpc_s(mpc_iter))
        end

        figure(F1);
        hmpc_s(mpc_iter) = plot(gam_s(1:gam_iter), max_setpoints_mpc(mpc_iter, 1:gam_iter),...
                            '-o', 'Color', colr, 'LineWidth', 2);

        upd(gam_iter)
        drawnow()

        
    end
 
    leg_str = sprintf('MPC, N=%.0f', N_mpc);
    set(hmpc_s(mpc_iter), 'DisplayName', leg_str);
    
    legend([hlin, hmpc_s]);
    xlabel('$\gamma$', 'interpreter', 'latex', 'FontSize', 16);
    ylabel('Max Ref', 'interpreter', 'latex', 'FontSize', 16);
    drawnow()
end
%%
% hmpc.DisplayName = 'MPC';
% legend([hlin, hmpc])

xlabel('$\gamma$', 'interpreter', 'latex', 'FontSize', 16)
ylabel('Max Ref', 'interpreter', 'latex', 'FontSize', 16)


ylim([0, 5.1])

max_ref_data.gam_s = gam_s;
max_ref_data.ref_range = ref_range;
max_ref_data.ref_step = ref_step;
max_ref_data.sim_struct = sim_struct;
max_ref_data.max_setpoints_mpc = max_setpoints_mpc;
max_ref_data.max_setpoints_lin = max_setpoints_linear;
title('with delay')
save('max_ref_data_dalay.mat', 'max_ref_data')
saveas(F1, 'max_ref_delay.fig')
