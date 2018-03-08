clc, clear
C_stage = 3.8e-6;
Imax = 100e-3;
Ts = 40e-6;


Vdiv = 2.56/(19.61+2.56);
% Vdiv_gain = 1/Vdiv; % from resistor measurement
Vdiv_gain = 10.6; % From 9v battery measurement
load( fullfile(PATHS.exp, 'x-axis_sines_info_out_2-8-2018-01.mat'))
Gpow_hat = modelFit.models.G_uz2pow;

Gpow = ss(Gpow_hat*Vdiv_gain);
Gpow.InputDelay = 0;
Gpow = minreal(Gpow*zpk([], [0 0], 1, Ts))


dVmax = (Ts/C_stage)*Imax
dVmax = 0.5;

G_stage = modelFit.models.G_pow2uz/Vdiv_gain;

G_stage.InputDelay = 0;
sys_old = modelFit.models.G_uz2stage;
sys = G_stage*Gpow;
% sys = modelFit.models.G_uz2stage;
F1 = figure(1);

[~,~,omegas] = bode(sys);
frfBode(Gpow, omegas/2/pi, F1,'-b', 'Hz');
frfBode(G_stage, omegas/2/pi, F1,'-k', 'Hz');
frfBode(sys_old, omegas/2/pi, F1,'-r', 'Hz');
frfBode(sys, omegas/2/pi, F1,'-g', 'Hz');

plotPZ_freqs(sys, F1)
plotPZ_freqs(G_stage, F1)
%%
figure(100)
pzplot(Gpow, 'b', G_stage, 'k')

figure(3)
pzplot(sys_old, 'g', sys, 'r')


%%
% sys_nodelay = sys;
sys_nodelay = modelFit.models.G_uz2stage;
sys_nodelay.InputDelay = 0;
% [sys_nodelay, gdrift] = eject_gdrift(sys_nodelay)
[wp_real_x, wz_real_x] = w_zp_real(sys_nodelay);
rho_1 = wz_real_x(1)/wp_real_x(1);

% zeta_x = [.9, .8, .6, .5 .5];
zeta_x = [.8, .7, .7, .5 .5, 0.25];
gams_x = [1.5, 1.5, 1.2, 1, 1, 1];
rhos_x = [rho_1*.99, 1, 1];
% rhos_x = [1, .5];

pint_x = 0.5*0;


pdes   = getCharDes(sys_nodelay, gams_x, pint_x, zeta_x, rhos_x, .25);
K_temp = place(sys_nodelay.a, sys_nodelay.b, pdes);
[Q0, R1, K_lqr] = inverseLQR(sys_nodelay, K_temp);
%     Q0 = blkdiag(Q0, zeros(Nd, Nd));

min(eig(Q0))
R1 = 10;

sys_recyc = SSTools.deltaUkSys(sys_nodelay);
Ns_mpc = size(sys_recyc.B, 1);
[Nx, Nu] = SSTools.getNxNu(sys_recyc);

Q1 = blkdiag(Q0, 0);
K_lqr = dlqr(sys_recyc.a, sys_recyc.b, Q1, R1);
Qp = dare(sys_recyc.a, sys_recyc.b, Q1, R1); 

sys_cl = ss(sys_recyc.a - sys_recyc.b*K_lqr, sys_recyc.b, sys_recyc.c, 0, sys.Ts);
figure(10)
pzplot(sys_recyc, sys_cl)
zgrid
xlim([0.7, 1])
ylim([-0.4, 0.4])
%%
% Now, try the state constraint with the condensed formulation.
x0_pow = [0;0];
du_max = 2.0;

R1 = 100; % OA works, QP fails
N_mpc = 100; % OA works, QP fails
N_mpc2 = 200;
ref = 4; %OA works, QP fails
x0 = -Nx*ref;

CON1 = CondenCon(Gpow, x0_pow, N_mpc);
CON1.add_state_con('box', dVmax);
% CON1.add_input_con('box', du_max);

CON1.xvec = x0_pow;

mpcProb1 = condensedMPCprob_OA(sys_recyc, N_mpc, Q1, Qp, R1);

mpcProb1.CON = CON1;

% Simulate the short MPC
sim_struct.K_lqr = K_lqr;
sim_struct.PLANT = sys_nodelay;
sim_struct.mpcProb1 = mpcProb1;
sim_struct.trun = 400*Ts;
sim_struct.du_max = du_max;
sim_struct.xss = Nx(1:end-1);
sim_struct.Nx = Nx;

sim_struct.uss_0 = 0;
sim_struct.Ts = Ts;
sim_struct.mpc_on = 1;
sim_struct.Gpow = Gpow;
sim_struct.x0_pow = x0_pow;


[Ympc, Umpc, dUmpc] = sim_MPC_fp(sim_struct, ref);
[ypow1, t1] = lsim(Gpow, dUmpc.Data, dUmpc.Time, x0_pow);


CON2 = CondenCon(Gpow, x0_pow, N_mpc2);
CON2.add_state_con('box', dVmax);
CON2.add_input_con('box', du_max);

mpcProb2 = condensedMPCprob_OA(sys_recyc, N_mpc2, Q1, Qp, R1);    
mpcProb2.CON = CON2; 
CON2.xvec = x0_pow;
[u2, X2] = mpcProb2.solve(x0, 'getX', 1);

t2 = [0:1:N_mpc2-1]'*sys.Ts;
ypow2 = lsim(Gpow, u2, t2, x0_pow);
ypow3 = Gpow.c*CON1.xvec;

figure(200);clf
    plot(ypow1)
    hold on
    plot(ypow2, '--')
    plot(ypow3, ':', 'LineWidth', 2)
    legend('MPC', 'CLQR')
    grid on
    title('Power Amplifer $\Delta y$')
    xlm = xlim;
    hold on;
    plot(xlm, [dVmax, dVmax], '--k')
    plot(xlm, -[dVmax, dVmax], '--k')
    
y2 = sys_recyc.c*(X2(:,1:end-1) - x0);


figure(9); clf
subplot(3,1,1)
    plot(Ympc.Time, Ympc.Data)
    hold on
    plot(t2, y2, '--')
    legend('MPC', 'CLQR')
    
    xlm = xlim;
    plot(xlm, [ref, ref]*1.01, ':k')
    plot(xlm, [ref, ref]*0.99, ':k')
    title('y(k)')
    grid on
subplot(3,1,2)
    hold on
    plot(dUmpc.Time, dUmpc.Data)
    plot(t2, u2, '--')
    grid on
    title('$\Delta u(k)$')

subplot(3,1,3)
    hold on
    plot(Umpc.Time, Umpc.Data)
    plot(t2, cumsum(u2), '--')
    grid on
    title('u(k)')


%%
% ----------------------------------------------------------------------- %
% -------------------- Now, try the time-optimal ------------------------ %
clc
k0 = 122
f_1 = CondensedTools.init_cond_resp_matrix(PHI_pow, 1, k0, C_pow);

H_1 = CondensedTools.zero_state_output_resp(Gpow, k0, 1);

binq = ones(k0,1)*dVmax - f_1*x0_pow*0;
sys_TO = eject_gdrift(sys_nodelay);

sys_TO = SSTools.deltaUkSys(sys_TO);
Nx_TO = SSTools.getNxNu(sys_TO);

x0_TO = Nx_TO*0;
xf = Nx_TO*ref;

addpath(genpath(fullfile(getMatPath, 'solvers/cvx')))
C=[];
du_max = 0.2;
for k=0:k0-1
    C = [sys_TO.a^k*sys_TO.b C];
end
e_k0 = zeros(1,k0);
e_k0(end) = 1; % Kth unit vector. [0 0 .... 0 0 1]. To pick out last element of u.
I = eye(length(sys_TO.b));

cvx_begin
    variable u(k0);
    variable t;
    L = sys_TO.a^k0*x0_TO + C*u;
    % minimize norm(L - xf)
    % minimize norm(u)
    minimize norm(t)
    subject to
%     norm(u, Inf) <= du_max;
    [H_1; -H_1]*u <= [binq; binq];
    % Enforce steady state requirement
    % xss = Axss + Buss --> (I-A)xss = Buss
    % (I-A)xss - Buss = 0. Not sure this is necessary.
    % If minimization succeeds, should be the case anyway.
    (I-sys_TO.a)*L - sys_TO.b*(e_k0*u) ==0;
    L - xf == 0;
cvx_end


results.U = u;
results.cvx_optval = cvx_optval;
results
%
u = [u; repmat(0, 400, 1)];

t = [0:1:length(u)-1]'*Ts;
y = lsim(sys_TO, u, t, x0_TO*0);
ypow = lsim(Gpow, u, t, x0_pow*0);
%
figure(9)
subplot(3,1,1), hold on
plot(t, y, '-.g')

subplot(3,1,2), hold on
plot(t, u, '-.g')
plot(t(k0), u(k0), 'x')
title('$\Delta u$')



subplot(3,1,3), hold on
plot(t, cumsum(u), '-.g')

title('u(k)')

figure(200)
hold on
plot(ypow)






%%
rmpath(genpath(fullfile(getMatPath, 'solvers/cvx')))






