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


G_stage = modelFit.models.G_pow2uz/Vdiv_gain;

G_stage.InputDelay = 0;
sys_old = modelFit.models.G_uz2stage;
sys = G_stage*Gpow;
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
sys_nodelay = sys;
sys_nodelay.InputDelay = 0;
% [sys_nodelay, gdrift] = eject_gdrift(sys_nodelay)
[wp_real_x, wz_real_x] = w_zp_real(sys);
rho_1 = wz_real_x(1)/wp_real_x(1);

% zeta_x = [.9, .8, .6, .5 .5];
zeta_x = [.8, .7, .7, .5 .5, 0.25];
gams_x = [1.5, 1.5, 1.2, 1, 1, 1];
% rhos_x = [rho_1*1.05, 1, 1];
rhos_x = [1, .5];
 
pint_x = 0.5*0;


%-----------------------------------------------------


P_x    = getCharDes(sys_nodelay, gams_x, pint_x, zeta_x, rhos_x, .25);
K_temp = place(sys_nodelay.a, sys_nodelay.b, P_x);
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

N_mpc = 300;
mpcProb0 = condensedMPCprob_OA(sys_recyc, N_mpc, Q1, Qp, R1);
% mpcProb0 = sparseMPCprob(sys_recyc, N_mpc, Q1, Qp, R1);
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
%     keyboard
    ydiffMat = [ydiffMat; ydiff_row];
end
%
ydiffMat = [zeros(size(ydiffMat,1), N_mpc), ydiffMat];

Ainq = [ydiffMat;
        -ydiffMat];
binq = ones(2*size(ydiffMat,1), 1)*dVmax;

% mpcProb0.Ainq = Ainq;
% mpcProb0.binq = binq;

[u, x] = mpcProb0.solve(x0, 'getX', 1);

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
figure(9);
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


