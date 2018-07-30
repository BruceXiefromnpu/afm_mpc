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
F1 = figure(1);

[~,~,omegas] = bode(sys);
frfBode(Gpow, omegas/2/pi, F1,'Hz', '-b');
frfBode(G_stage, omegas/2/pi, F1,'Hz', '-k');
frfBode(sys_old, omegas/2/pi, F1,'Hz', '-r');
frfBode(sys, omegas/2/pi, F1,'Hz', '-g');

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
[wp_real_x, wz_real_x] = w_zp_real(sys);
rho_1 = wz_real_x(1)/wp_real_x(1);

% zeta_x = [.9, .8, .6, .5 .5];
zeta_x = [.8, .7, .7, .5 .5, 0.25];
gams_x = [1.5, 1.5, 1.2, 1, 1, 1];
% rhos_x = [rho_1*1.05, 1, 1];
rhos_x = [1, 1];

pint_x = 0.5*0;


pdes   = getCharDes(sys_nodelay, gams_x, pint_x, zeta_x, rhos_x, .25);
K_temp = place(sys_nodelay.a, sys_nodelay.b, pdes);
[Q0, R1, S0, K_lqr] = inverseLQR_cross_weight(sys_nodelay, K_temp);
%     Q0 = blkdiag(Q0, zeros(Nd, Nd));

min(eig(Q0))
R1 = 10;

pls1 = sort(eig(sys_nodelay.a - sys_nodelay.b*K_temp))
pls2 = sort(eig(sys_nodelay.a - sys_nodelay.b*K_lqr))
%%
[Q1, R1, S1, K_lqr2] = inverseLQR_cross_weight(sys_nodelay, K_lqr);


K_lqr
K_lqr2
K_temp


%%
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

N_mpc = 300;

mpcProb0 = sparseMPCprob(sys_recyc, N_mpc, Q1, Qp, R1);
ref = 5;
x0 = -Nx*ref;
% [u, x] = mpcProb0.solve(x0, 'getX', 1);

Cpowpick = [zeros(1, 12), Gpow.c, 0];
% now, we should be able to construct a diff matrix that will pick out the 
Ns = 15;
ydiffMat = []
for kk=1:N_mpc-1
    ydiff_row = [zeros(1, (kk-1)*Ns), -Cpowpick, Cpowpick, zeros(1, N_mpc*Ns - kk*Ns)];
%     ydiff_row = [zeros(1, (kk-1)*Ns), Cpowpick, zeros(1, N_mpc*Ns+Ns - kk*Ns)];
    ydiffMat = [ydiffMat; ydiff_row];
end

ydiffMat = [zeros(size(ydiffMat,1), N_mpc), ydiffMat];

Ainq = [ydiffMat;
        -ydiffMat];
binq = ones(2*size(ydiffMat,1), 1)*dVmax;

mpcProb0.Ainq = Ainq;
mpcProb0.binq = binq;

[u, x] = mpcProb0.solve(x0, 'getX', 1);
%%
clc
figure(200); clf
    xpow = x(13:14,:);
    ypow = Gpow.c*xpow;
    plot(diff(ypow))
    grid on
    title('Power Amplifer $\Delta y$')
    xlm = xlim;
    hold on;
    plot(xlm, [dVmax, dVmax], '--k')
    plot(xlm, -[dVmax, dVmax], '--k')

y = sys_recyc.c*(x(:,1:end-1) - x0);

t = [0:1:N_mpc-1]'*sys.Ts;
figure(9); clf
subplot(3,1,1)
    plot(t, y)
    hold on
    xlm = xlim
    plot(xlm, [ref, ref]*1.01, ':k')
    plot(xlm, [ref, ref]*0.99, ':k')

grid on
subplot(3,1,2)
plot(t, u)
grid on

subplot(3,1,3)
plot(t, cumsum(u))
grid on

%%
% Now, try the state constraint with the condensed formulation.
clc
N_mpc = length(u);
% N_mpc = 4;
% Y = f * x0 + H*U, where f,H are matrices, and U is a stacked vector 
% [u0, .... u_{N-1}]'
U = u(1:N_mpc);
U = U(:);
if 0
    Gpow_recyc = SSTools.deltaUkSys(Gpow);
    x0_pow = x0(13:15);
    PHI_pow = Gpow_recyc.a;
    C_pow = Gpow_recyc.c;

    Ns = 3; nu = 1;
    no=1;

    f_0 = CondensedTools.init_cond_resp_matrix(PHI_pow, 0, N_mpc-no, C_pow);
    f_1 = CondensedTools.init_cond_resp_matrix(PHI_pow, 1, N_mpc, C_pow);

    H_0 = CondensedTools.zero_state_output_resp(Gpow_recyc, N_mpc-1, 0);

    H_1 = CondensedTools.zero_state_output_resp(Gpow_recyc, N_mpc, 1);
    % We have dY = Y1 - Y0
    %            = f1*x0 + H1*U - [f0*x0 + H0*U]
    %            = (f1 - f0)*x0 + (H1-H0)*U 
    %            <= dyMAX
    df = f_1 - f_0;
    dH = H_1 - H_0;

    binq = ones(N_mpc,1)*dVmax;
    binq = binq - df*x0_pow;
    Ainq = dH;
    y1pow_fromMat = [f_1*x0_pow + H_1*( U)];
    y0pow_fromMat = [f_0*x0_pow + H_0*( U)];
    dy_fromMat = y1pow_fromMat - y0pow_fromMat;
else
    % x0_pow = x0(13:14);
    x0_pow = [0;0];
    PHI_pow = Gpow.a;
    C_pow = Gpow.c;

    Ns = 2; nu = 1;
    no=1;

    f_0 = CondensedTools.init_cond_resp_matrix(PHI_pow, 1, N_mpc, C_pow);
    H_0 = CondensedTools.zero_state_output_resp(Gpow, N_mpc, 1);

    binq = ones(N_mpc,1)*dVmax - f_0*x0_pow*0;
    Ainq = H_0;
    dy_fromMat = [f_0*x0_pow*0 + H_0*( U)];
end

     
%%
figure(300), clf
hold on
plot(diff(ypow))
plot(dy_fromMat, '--')
% plot(y1pow_fromMat, '-.')
    
mpcProb2 = condensedMPCprob_OA(sys_recyc, N_mpc, Q1, Qp, R1);    


mpcProb2.Ainq = Ainq;
mpcProb2.lbAinq = -binq;
mpcProb2.ubAinq = binq;

[u2, X2] = mpcProb2.solve(x0, 'getX', 1);

%%
figure(200); 
    xpow2 = X2(13:14,:);
    ypow2 = Gpow.c*xpow2;
    plot(diff(ypow2), '--')
    grid on
    title('Power Amplifer $\Delta y$')
    xlm = xlim;
    hold on;
    plot(xlm, [dVmax, dVmax], '--k')
    plot(xlm, -[dVmax, dVmax], '--k')

y2 = sys_recyc.c*(X2(:,1:end-1) - x0);

t = [0:1:N_mpc-1]'*sys.Ts;
figure(9);
subplot(3,1,1)
    plot(t, y2, '--')
    hold on
    xlm = xlim
    plot(xlm, [ref, ref]*1.01, ':k')
    plot(xlm, [ref, ref]*0.99, ':k')

grid on
subplot(3,1,2)
hold on
plot(t, u2, '--')
grid on

subplot(3,1,3)
hold on
plot(t, cumsum(u2), '--')
grid on


%%
% ----------------------------------------------------------------------- %
% -------------------- Now, try the time-optimal ------------------------ %
clc
k0 = 150
f_0 = CondensedTools.init_cond_resp_matrix(PHI_pow, 0, k0-no, C_pow);
f_1 = CondensedTools.init_cond_resp_matrix(PHI_pow, 1, k0, C_pow);
df = f_1 - f_0;

H_0 = CondensedTools.zero_state_output_resp(Gpow_recyc, k0-1, 0);
H_1 = CondensedTools.zero_state_output_resp(Gpow_recyc, k0, 1);
dH = H_1;
% - H_0;

binq = ones(k0,1)*dVmax - df*x0_pow*0;

xf = x0*0;
addpath(genpath(fullfile(getMatPath, 'solvers/cvx')))
C=[];
du_max = 0.2;
for k=0:k0-1
    C = [sys_recyc.a^k*sys_recyc.b C];
end
e_k0 = zeros(1,k0);
e_k0(end) = 1; % Kth unit vector. [0 0 .... 0 0 1]. To pick out last element of u.
I = eye(length(sys_recyc.b));

cvx_begin
    variable u(k0);
    variable t;
    L = sys_recyc.a^k0*x0 + C*u;
    % minimize norm(L - xf)
    % minimize norm(u)
    minimize norm(t)
    subject to
%     norm(u, Inf) <= du_max;
    [dH; -dH]*u <= [binq; binq];
    % Enforce steady state requirement
    % xss = Axss + Buss --> (I-A)xss = Buss
    % (I-A)xss - Buss = 0. Not sure this is necessary.
    % If minimization succeeds, should be the case anyway.
    (I-sys_recyc.a)*L - sys_recyc.b*(e_k0*u) ==0;
    L - xf == 0;
cvx_end


results.U = u;
results.cvx_optval = cvx_optval;
results
%%
u = [u; repmat(0, 400, 1)];

t = [0:1:length(u)-1]'*Ts;
y = lsim(sys_recyc, u, t, x0*0);
ypow = lsim(Gpow_recyc, u, t, x0_pow*0);
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
plot(diff(ypow))






%%
rmpath(genpath(fullfile(getMatPath, 'solvers/cvx')))






