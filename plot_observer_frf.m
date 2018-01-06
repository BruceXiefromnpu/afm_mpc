clear;
clc;

addpath('functions')
% ---------- Load Parametric Models  -----------
modFitPath    = 'x-axis_sines_info_out_12-10-2017-04.mat';

load(fullfile(PATH_sysid, modFitPath));



sys = modelFit.sys_x_eigen;


if 1
    Nd =sys.InputDelay;
    sys = absorbDelay(sys)
else
    sys.InputDelay = 0
    Nd = 0;
end
   

Ts = sys.Ts;
if 0
    w_s = modelFit.frf.w_s;
    sys_frf = modelFit.frf.G_frf;
else
   [sys_frf, w_s] = getFRF(sys); 
end
freqs = w_s/2/pi;


% -----------------------------------------------------------------------%
% Design K:
pstyle = 'b';
gam_s = [1., 1., 1, 1, 1]; % factors to increase cmplx mode freqs by
alp_s = [.7 .7 .4 .4, .2]; % desired damping of cmplx modes
rho_s = [1.2, 1]; % factors to shift real modes by


% fdbk is a class which computes the gain and also stores the data used to
% generate K with. 
Kfdbk = fdbk(sys, 'gams', gam_s, 'alps',...
         alp_s, 'rhos', rho_s, 'doDelay', 1, 'rad', 0.1 );
K = Kfdbk.K;      

[Q, R, K] = inverseLQR(sys, K);
% Q = 10*sys.c'*sys.c;
K = dlqr(sys.a, sys.b, Q, R);
Nbar = getNbar(sys, K);
syscl = ss(sys.a - sys.b*K, sys.b, sys.c, 0, Ts);
figure(20); clf;
pzplot(syscl)

% -----------------------------------------------------------------------%
%                        Design L                                        %
% -----------------------------------------------------------------------%

BvB = 1*modelFit.sys_x_eigen.b(1:12)*modelFit.sys_x_eigen.b(1:12)';
Qdelay = zeros(Nd, Nd)+0.0001;
% Qdelay = sys.b(13:end)*sys.b(13:end)'*.00000;
Qw = blkdiag(modelFit.Qw + BvB, Qdelay);

% Qw = sys.c'*sys.c*10000;

Rw = modelFit.Rw;
L = dlqr(sys.a', sys.c', Qw, Rw)';
syscl_obs = ss(sys.a - L*sys.c, sys.b, sys.c, 0, Ts);
figure(30);clf
pzplot(syscl_obs)
hold on
pzplot(sys, 'r')

% Create ss systems.
A_tilde = sys.a - sys.b*K - L*sys.c;
D1 = ss(A_tilde, L, K, 0, Ts);
D2 = ss(A_tilde, sys.b, K, 0, Ts);


% Pure state feedback loop gain:
LL = ss(sys.a, sys.b, K, 0, Ts);

% Estimated state feedback loop gain:
LL_obs = D1*sys;

H = Nbar*sys*(1-D2)/(1 + sys*D1);

D1_frf = getFRF(D1, w_s);
D2_frf = getFRF(D2, w_s);
H_frf = Nbar*sys_frf.*(1 - D2_frf)./(1 + sys_frf.*D1_frf);
LL_frf = getFRF(LL, w_s);
LL_obs_frf = D1_frf.*sys_frf;

% H_frf = getFRF(H, w_s);



h = [];
% gobjects(4,1)
F1 = figure(1); clf
h_i = frfBode(H_frf, freqs, F1, '-k', 'Hz');
h_i.DisplayName = 'H closed loop';
h = [h; h_i];

% h(2) = frfBode(sys*D1, freqs, F1, '--k', 'Hz');

h_i = frfBode(LL_frf, freqs, F1, '--r', 'Hz');
h_i.DisplayName = 'L = K(zI-A)^-1B';
h = [h; h_i];

h_i = frfBode(LL_obs_frf, freqs, F1, '--b', 'Hz');
h_i.DisplayName = 'L = G(s)K(zI - A + BK + LC)L';
h = [h;h_i];
% h(2).DisplayName = 'G*D1';

legend(h, 'location', 'northwest')

figure(20)
hold on;
pzplot(LL, 'r')
[x, kk]=rlocus(LL);
for k=1:size(x,1)
    plot_ind(real(x(k, :)), imag(x(k, :)))
end
    
% %%
% g = zpk([-1+20j, -1-20j], [-5, -1+10j, -1-10j, -2+20j, -2-20j], 1)
% 
% g = c2d(g, .001);
% 
% sys = ss(g);
% K = dlqr(sys.a, sys.b, eye(5), 10);
% figure(2)
% bode(sys)

if 0
    zrs = tzero(LL)
    kk = find(abs(zrs)>1)

    zrs(kk) = 0.8;

    [C, D] = place_zeros(sys, zrs');
    K2 = dlqr(sys.a, sys.b, C'*C, 1, S);
    LL_placed = ss(sys.a, sys.b, K2, D, Ts);


    figure(30); clf
    hold on;
    pzplot(syscl, 'b')
    pzplot(LL_placed, 'r')
    [x, kk]=rlocus(LL_placed);
    for k=1:size(x,1)
        plot_ind(real(x(k, :)), imag(x(k, :)))
    end

end

% Nyquist for pure state feedback loop gain
figure(1000); clf;

nyquist_frf(LL_frf, w_s)
hold on
t = 0.001:0.01:2*pi;
% The Nyquist diagram of discrete LQR Loop is guaranteed to be outside of
% the circle centered at (-1, 0j) with radius 1/(B'*PB + 1)^1/2, P is 
% solution to DARE. THis is in contrast to being outside of a circle of 
% radius one as in continuous time.

P = dare(LL.a, LL.b, Q, R);

rad = 1/(sys.b'*P*sys.b + 1)^(0.5);

x = rad*sin(t);
y = rad*cos(t);

plot(x-1,y, '--')
plot(x-1,y, '--')
xlim([-1.5, -0.1])
ylim([-1., 1])
grid
title('Nyquist of Loop gain pure state feedback')

% Nyquist for ESTIMATE state feedback loop gain
figure(1001); clf;

nyquist_frf(LL_obs_frf, w_s)
hold on
t = 0.001:0.01:2*pi;
% The Nyquist diagram of discrete LQR Loop is guaranteed to be outside of
% the circle centered at (-1, 0j) with radius 1/(B'*PB + 1)^1/2, P is 
% solution to DARE. THis is in contrast to being outside of a circle of 
% radius one as in continuous time.

P = dare(LL.a, LL.b, Q, R);

rad = 1/(sys.b'*P*sys.b + 1)^(0.5);

x = rad*sin(t);
y = rad*cos(t);

plot(x-1,y, '--')

xlim([-1.5, -0.1])
ylim([-1., 1])
grid on
title('Nyquist of Loop gain with observer')
