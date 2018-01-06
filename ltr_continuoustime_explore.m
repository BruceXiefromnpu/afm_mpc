
clear;
% close all
clc
figpath = '/home/arnold/rabraker.com/content/software/figures/';
m = 2;
b = 2;
bo = 1/m
k = 10;

w = sqrt(k/m);
gamma = b/m;


G = tf(bo, [1 gamma, w^2])

wp2 = 7.8;
wz2 = 7.;
G2 = tf([1 2*0.02*wz2 wz2^2], [1 2*0.031*wp2, wp2^2])*(wp2^2/wz2^2);
G = zpk(G*G2)*10;
sys = ss(G);

[sys_frf, w_s] = getFRF(sys);


% figure(2);clf; hold on
% qs = logspace(log(1), log10(2^4), 4);
% for k = 1:length(qs)
% q = qs(k);
% Q = q*sys.c'*sys.c;
% K = lqr(sys.a, sys.b, Q, 1);
% syscl=ss(sys.a-sys.b*K, sys.b, sys.c,0);
% pzplot(syscl, 'k')
% drawnow()
% 
% end

Q = diag([1, 100, 0 0.0001]);
% Q = diag([1, 100]);
K = lqr(sys.a, sys.b, Q, 1);


L_pure = ss(sys.a, sys.b, K, 0);
L_pure_frf = getFRF(L_pure, w_s);

Qw = sys.b*sys.b'*10;

K_e = lqr(sys.a', sys.c', Qw, 1)';

A_tilde = sys.a - sys.B*K - K_e*sys.c;

D1 = ss(A_tilde, K_e, K, 0);
D2 = 1-ss(A_tilde, sys.b, K, 0);
D2_frf = getFRF(D2, w_s);

Nbar = getNbar(sys, K);

H = Nbar*sys*D2/(1 + G*D1);

H_frf = getFRF(H, w_s);

L_obs = sys*D1;
L_obs_frf = getFRF(L_obs, w_s);

figure(100);clf
rlocus(L_obs)

% ----------------------------------------------------------------------- %
%                       Bode Plots                                        %
freqs = w_s/2/pi;
F1 = figure(200);clf
h1 = frfBode(sys_frf, freqs, F1, '-k', 'Hz');
h1.DisplayName = 'G(s)';
h2 = frfBode(L_pure_frf, freqs, F1, 'b', 'Hz');
h2.DisplayName = 'K(sI - A)^-1B';
h3 = frfBode(L_obs_frf, freqs, F1, '--r', 'Hz');
h3.DisplayName = 'K(sI - A + BK +LC)^-1BG(s)';
h4 = frfBode(H_frf, freqs, F1, ':k', 'Hz');
h4.DisplayName = 'Closed Loop (obs)';
legend([h1, h2, h3, h4], 'location', 'southwest')

% ----------------------------------------------------------------------- %
%                       Nyquist for state feedback loop gain              %
F1000 = figure(1000); clf
nyquist_frf(L_pure_frf, w_s)
hold on
t = 0.001:0.01:2*pi;
x = sin(t);
y = cos(t);

plot(x-1,y, '--')
plot(x-1,y, '--')
% xlim([-1.5, -0.1])
% ylim([-1., 1])
grid
title('Nyquist of Loop gain pure state feedback')
saveon = 0
if saveon
    saveas(F1000, fullfile(figpath, 'lqr_nyquist_stateFB_example01.svg'))
end
    
% ----------------------------------------------------------------------- %
%             Nyquist for estimator state feedback loop gain              %
F1001=figure(1001); clf
nyquist_frf(L_obs_frf, w_s)
hold on
t = 0.001:0.01:2*pi;
x = sin(t);
y = cos(t);

plot(x-1,y, '--')
plot(x-1,y, '--')
% xlim([-1.5, -0.1])
% ylim([-1., 1])
grid
title('Nyquist of Loop gain estimator feedback')


if saveon
    saveas(F1001, fullfile(figpath, 'figures/lqr_nyquist_obsFB_example01.svg'))
end

%% ---------------------------------------------------------------------- %
% Now, let's try this in discrete time and see if the LTR works the same
% way.
clc
Ts = 0.05;
sys_z = c2d(sys, Ts);


Qz = sys_z.c'*sys_z.c*100;
% Qz = diag([10.0, 200, 10 0.01]);
% Qz = diag([100.0, 2000]);
% Qz = rand(4,4); Qz = (Qz+Qz')/2;
% Qz = Qz + 2*abs(min(eig(Qz)))*eye(4);
Rz = 1;
[Kz, Pz] = dlqr(sys_z.a, sys_z.b, Qz, Rz);
Kz=Kz;

Lz_pure = ss(sys_z.a, sys_z.b, Kz, 0, Ts);
[Lz_pure_frf,w_s] = getFRF(Lz_pure);
% [sys_z_frf, w_s] = getFRF(sys_z);


F2000 = figure(2000); clf
nyquist_frf(Lz_pure_frf, w_s)
hold on

t = 0.001:0.01:2*pi;
x = sin(t);
y = cos(t);

gamma = 1/(sys_z.b'*Pz*sys_z.b + 1)^(0.5);
rad = gamma;


plot(rad*x-1, rad*y)
plot(x-1, y, '--')
grid
%%
% ----------------------------------------------------------------------- %
% Now design estimator.

Qwz = sys_z.b*sys_z.b'*100;

Kz_e = dlqr(sys_z.a', sys_z.c', Qwz, 1)';

Az_tilde = sys_z.a - sys_z.B*Kz - Kz_e*sys_z.c;

D1z = ss(Az_tilde, Kz_e, Kz, 0, Ts);
D2z = 1-ss(Az_tilde, sys_z.b, Kz, 0, Ts);
D2z_frf = getFRF(D2z, w_s);

Nbar = getNbar(sys_z, Kz);

Hz = Nbar*sys_z*D2z/(1 + sys_z*D1z);

Hz_frf = getFRF(Hz, w_s);

Lz_obs = sys_z*D1z;
Lz_obs_frf = getFRF(Lz_obs, w_s);



F3000 = figure(3000); clf

nyquist_frf(Lz_obs_frf, w_s)
hold on

t = 0.001:0.01:2*pi;
x = sin(t);
y = cos(t);



plot(rad*x-1, rad*y)
plot(x-1, y, '--')
grid on
% xlim([-2, 2])
% ylim([-2, 2])

% ----------------------------------------------------------------------- %
%                       Bode Plots                                        %
% freqs = w_s/2/pi;
% F10 = figure(10);clf
% h1 = frfBode(sys_z_frf, freqs, F10, '-k', 'Hz');
% h1.DisplayName = 'G(s)';
% h2 = frfBode(Lz_pure_frf, freqs, F10, 'b', 'Hz');
% h2.DisplayName = 'K(sI - A)^-1B';
% h3 = frfBode(Lz_obs_frf, freqs, F10, '--r', 'Hz');
% h3.DisplayName = 'K(sI - A + BK +LC)^-1BG(s)';
% h4 = frfBode(Hz_frf, freqs, F10, ':k', 'Hz');
% h4.DisplayName = 'Closed Loop (obs)';
% legend([h1, h2, h3, h4], 'location', 'southwest')

% ----------------------------------------------------------------------- %
%                       Nyquist for state feedback loop gain              %

%%


% Pz = dare(Lz_pure.a, Lz_pure.b, Qz, Rz);
% gamma = 1/sqrt(sys_z.b'*Pz*sys_z.b + 1)
% 
% plot(x-1,y, '--')
% 
% plot(gamma*x-1, gamma*y, ':r')

% xlim([-1.5, -0.1])
% ylim([-1., 1])
grid
title('Discrete Nyquist of Loop gain pure state feedback')
saveon = 0;
if saveon
    saveas(F2000, fullfile(figpath, 'lqr_nyquist_stateFB_example01.svg'))
end
%%    
% ----------------------------------------------------------------------- %
%             Nyquist for estimator state feedback loop gain              %
F2001=figure(2001); clf
nyquist_frf(Lz_obs_frf, w_s)
hold on
t = 0.001:0.01:2*pi;
x = sin(t);
y = cos(t);

plot(x-1,y, '--')
plot(x-1,y, '--')
% xlim([-1.5, -0.1])
% ylim([-1., 1])
grid
title('Diszrete Nyquist of Loop gain estimator feedback')


if saveon
    saveas(F2001, fullfile(figpath, 'figures/lqr_nyquist_obsFB_example01.svg'))
end



