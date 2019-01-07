clear;
clc;

addpath('functions')
% ---------- Load Parametric Models  -----------
modFitPath    = 'x-axis_sines_info_out_12-10-2017-04.mat';

Ts = 0.01;
sys = c2d(ss(tf(100, [1, 20, 100])), Ts);

Q = eye(2);
R = 1;
K = dlqr(sys.a, sys.b, Q, R);
L = dlqr(sys.a', sys.c', Q, R)';


% Create ss systems.
A_tilde = sys.a - sys.b*K - L*sys.c;
D1 = ss(A_tilde, L, K, 0, Ts);
D2 = ss(A_tilde, sys.b, K, 0, Ts);
Nbar = SSTools.getNbar(sys, K);

% Pure state feedback loop gain:
LL = ss(sys.a, sys.b, K, 0, Ts);

% Estimated state feedback loop gain:
LL_obs = D1*sys;

H = Nbar*sys*(1-D2)/(1 + LL_obs);

S = feedback(1, LL_obs);

h = [];
F1 = figure(1); clf
% h_i = bode(H)
% h_i.DisplayName = 'H closed loop';
bode(H)
hold on
bode(LL, '--r')
%h_i.DisplayName = '$L = K(zI-A)^{-1}B$';

bode(LL_obs, '--b')
bode(S)

% h_i.DisplayName = '$L = G(s)K(zI - A + BK + LC)L$';

legend('H-CL', '$L = K(zI-A)^{-1}B$', ['$L = G(s)K(zI - A + BK + LC)L$'], 'S');


figure(2); clf
step(S)
hold on
step(H)
step(sys*S)
% %%
% % Nyquist for pure state feedback loop gain
% figure(1000); clf;

% nyquist_frf(LL_frf, w_s)
% hold on
% t = 0.001:0.01:2*pi;
% % The Nyquist diagram of discrete LQR Loop is guaranteed to be outside of
% % the circle centered at (-1, 0j) with radius 1/(B'*PB + 1)^1/2, P is 
% % solution to DARE. THis is in contrast to being outside of a circle of 
% % radius one as in continuous time.

% P = dare(LL.a, LL.b, Q, R);

% rad = sqrt(R)/sqrt(sys.b'*P*sys.b + R);

% x = rad*sin(t);
% y = rad*cos(t);

% plot(x-1,y, '--')
% plot(x-1,y, '--')
% % xlim([-1.5, -0.1])
% % ylim([-1., 1])
% grid
% title('Nyquist of Loop gain pure state feedback')

% % Nyquist for ESTIMATE state feedback loop gain
% figure(1001); clf;

% nyquist_frf(LL_obs_frf, w_s)
% hold on
% t = 0.001:0.01:2*pi;
% % The Nyquist diagram of discrete LQR Loop is guaranteed to be outside of
% % the circle centered at (-1, 0j) with radius 1/(B'*PB + 1)^1/2, P is 
% % solution to DARE. THis is in contrast to being outside of a circle of 
% % radius one as in continuous time.

% P = dare(LL.a, LL.b, Q, R);

% rad = 1/(sys.b'*P*sys.b + 1)^(0.5);

% x = rad*sin(t);
% y = rad*cos(t);

% plot(x-1,y, '--')

% % xlim([-1.5, -0.1])
% % ylim([-1., 1])
% grid on
% title('Nyquist of Loop gain with observer')
