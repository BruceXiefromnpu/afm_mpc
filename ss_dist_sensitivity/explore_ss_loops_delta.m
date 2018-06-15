clear;
clc;

addpath('functions')
% ---------- Load Parametric Models  -----------
modFitPath    = 'x-axis_sines_info_out_12-10-2017-04.mat';

Ts = 0.01;
sys = c2d(ss(tf(100, [1, 20, 100])), Ts);

Q = eye(3);
R = 1;
sys_recyc = SSTools.deltaUkSys(sys);

K_lqr = dlqr(sys_recyc.a, sys_recyc.b, Q, R);
Kx = K_lqr(1:end-1);
Nbar = SSTools.getNbar(sys_recyc, K_lqr);
% Nbar = 1;

Nx = SSTools.getNxNu(sys_recyc);

Ku = K_lqr(end);

L = dlqr(sys.a', sys.c', sys.b*sys.b', R)';

% Create ss systems.
A_tilde = sys.a - sys.b*Kx - L*sys.c;

AA_ = [A_tilde, sys.b-sys.b*Ku;
       -Kx,     1-Ku];
BB_ = [sys.b*Nbar; Nbar];
LL_ = [L; 0];


K_c = K_lqr;
K_c(end) = K_c(end)-1;
D1 = ss(AA_, LL_, K_c, 0, Ts);
D2 = ss(AA_, BB_, -K_c, Nbar, Ts);



% Estimated state feedback loop gain:
loop = minreal(sys*D1);

H_yr = (sys*D2/( 1 + loop));
D3 = ss(AA_, LL_+BB_, K_c, -Nbar, Ts);
H_err = (1 + sys*D3)/(1 + loop);



PLANT = sys;
sys_obs = sys;
ref = 1;
dist = 0;
x0 = PLANT.b*0;
x0_obs = sys_obs.b*0;
trun = 2;
sim('AFMss_obshas_uk')

u = Y.Time*0+ref;
y_H = lsim(H_yr, u, Y.Time);
e_H = lsim(H_err, u, Y.Time);


figure(2); clf
plot(Y.Time, Y.Data)
% step(S)
hold on
plot(Y.Time, y_H, '--')
% step(sys*S)


figure(3); clf
plot(Y.time, ref-Y.data)
hold on
plot(Y.Time, e_H, '--')


pole(H_yr)


loop = minreal(sys*D1);
M = [loop.a, loop.b;
     loop.c, 1];

P = eye(size(M,1));
P(end,end) = 0;
pp = eig(M, P)
