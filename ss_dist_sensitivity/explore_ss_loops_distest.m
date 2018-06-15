clear;
clc;

addpath('functions')
% ---------- Load Parametric Models  -----------
modFitPath    = 'x-axis_sines_info_out_12-10-2017-04.mat';

Ts = 0.01;
sys = c2d(ss(tf(100, [1, 2, 100])), Ts);

% Q = 1100*eye(3);
% R = 1;

Px = getCharDes(sys, [1, 1], 0.9, [0.7], [])
[chat, dhat] = place_zeros(sys, Px);
Q = chat'*chat;
S = chat'*dhat;
R = dhat^2+0;

K_lqr = dlqr(sys.a, sys.b, Q, R, S);

Kx = K_lqr;
Nbar = SSTools.getNbar(sys, K_lqr);
% Nbar = 1;

Nx = SSTools.getNxNu(sys);

Ku = K_lqr(end);

Qw = sys.b*sys.b'*10;
Lx = dlqr(sys.a', sys.c', Qw, 1)';
p_int_d = 0.8;
[L, sys_obs, Ident_obs, C_ydist] = DistEst.output_dist_est(sys, Lx, p_int_d);
Ld = L(end);
Lx = L(1:end-1);
Cd = sys_obs.c(end)

% Create ss systems.
A_tilde = sys.a - sys.b*Kx - Lx*sys.c;

AA_ = [A_tilde,  -(Nbar*sys.b+Lx*Cd);
      -Ld*sys.c,        1-Ld*Cd];

BB_ = [sys.b*Nbar; 0];
LL_ = [Lx; Ld];


K_c = [K_lqr, Nbar];
D1 = ss(AA_, LL_, K_c, 0, Ts);
D2 = ss(AA_, BB_, -K_c, Nbar, Ts);



% % Estimated state feedback loop gain:
loop = sys*D1;

S = 1/(1 + loop);
H_yr = sys*D2*S;

H_yd = sys*S;

F1 = figure(1); clf
[~,~, omegas] = bode(H_yr);

h_i = frfBode(H_yr, omegas, F1, '-b', 'rad');
hold on
h_S = frfBode(S, omegas, F1, '--k', 'rad');

h_i.DisplayName = 'H closed loop';
h_S.DisplayName = 'Sensitivity';
legend([h_i, h_S])
% bode(LL, '--r')
% %h_i.DisplayName = '$L = K(zI-A)^{-1}B$';

% bode(LL_obs, '--b')
% bode(S)

% % h_i.DisplayName = '$L = G(s)K(zI - A + BK + LC)L$';

% legend('H-CL', '$L = K(zI-A)^{-1}B$', ['$L = G(s)K(zI - A + BK + LC)L$'], 'S');

PLANT = sys;

ref = 0;
dist = 1;
x0 = PLANT.b*0;
x0_obs = sys_obs.b*0;
trun = 2;
sim('AFMss_dist_obshas_uk')

u = Y.Time*0+ref;
d = 0*u + dist;
y_H = lsim(H_yr, u, Y.Time);
y_Hyd = lsim(H_yd, d, Y.Time);
% e_H = lsim(H_err, u, Y.Time);


figure(2); clf
plot(Y.Time, Y.Data)
% step(S)
hold on
% plot(Y.Time, y_H, '--')
plot(Y.Time, y_Hyd, '--')

return
figure(3); clf
plot(Y.time, ref-Y.data)
hold on
plot(Y.Time, e_H, '--')


pole(H)


loop = minreal(sys*D1);
M = [loop.a, loop.b;
     loop.c, 1];

P = eye(size(M,1));
P(end,end) = 0;
pp = eig(M, P)
