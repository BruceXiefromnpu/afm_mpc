
addpath(fullfile(getMatPath(), 'afm_mpc_journal', 'functions'))
addpath(fullfile(getMatPath(), 'afm_mpc_journal', 'functions', 'canon'))
addpath('~/gradschool/sysID/matlab/functions/')

clear, clc
saveon = true;
plants = CanonPlants.plants_ns14(9,2);
G = plants.SYS;
G_recyc = plants.sys_recyc;
Ts = G_recyc.Ts;

pmgm_figfile = 'PMGM_vs_gamma_both.svg';
sens_ts_figfile_both = 'GainS_TS_vs_gamma_both.svg';
bw_ts_figfile_both = 'BW_TS_vs_gamma_both.svg';

Qw = plants.SYS.b*plants.SYS.b'*50;
Lx = G.a*dlqr(plants.SYS.a', plants.SYS.c', Qw, 1)';
p_int_d = 0.7;
[LxLd, G_obsDist, Ident_obs, C_ydist] = DistEst.output_dist_est(G, Lx, p_int_d);

cmplx_rad = 0.9;
% Constant sigma LQR weights
[Q_cs, R_cs, S_cs, P_x] = build_control_constsigma(G_recyc, cmplx_rad);
% Chooze zeta LQR weights
can_cntrl = CanonCntrlParams_ns14();
[Q_cz, R_cz, S_cz] = build_control(G_recyc, can_cntrl);

gam = 25;
Ts = G.Ts;
omegas = logspace(log10(1), log10(12500*2*pi), 1000);

K_cz = dlqr(G_recyc.a, G_recyc.b, Q_cz, R_cz+gam, S_cz);
[S_cz, ~, Hyr_cz, ~, L_cz] = ss_loops_delta_dist(G, G_recyc, G_obsDist, K_cz, LxLd);

K_cs = dlqr(G_recyc.a, G_recyc.b, Q_cs, R_cs+gam, S_cs);
[S_cs, ~, Hyr_cs, ~, L_cs] = ss_loops_delta_dist(G, G_recyc, G_obsDist, K_cs, LxLd);

[gm_cz, pm_cz] = margin(L_cz);
[gm_cs, pm_cs] = margin(L_cs);



figure(1)
nyquist(L_cz)

figure(2)
nyquist(L_cs)
%%
[re, im, w] = nyquist(L_cz);
[gm, pm, wcg, wcp] = margin(L_cz)
gm_abs = 10^(gm/20)
[re, im, w] = nyquist(L_cz, sort([w; wcg; wcp]));

idx_wcg = w(w==wcg);
idx_wcp = w(w==wcp);
L_wcp = squeeze(freqresp(L_cz, wcp));
L_wcg = squeeze(freqresp(L_cz, wcg));


re = re(:);
im = im(:);

t = [0:0.01:2*pi]';
x = cos(t);
y = sin(t);

figure(3); clf
subplot(2,2,1)
plot(re, im)
hold on
grid on
plot(re, -im, '--')
plot(x,y, ':k')
legend('$\omega >0$', '$\omega<0$')
title('Nyquist diagriam')

subplot(2,2,2)
plot(re, im)
hold on
grid on
plot(re, -im, '--')
plot(x,y, ':k')
plot([-1], 0, 'rx', 'MarkerSize', 5)
plot([0 real(L_wcp)], [0, imag(L_wcp)], '-k')
plot([real(L_wcp)], [imag(L_wcp)], '.', 'MarkerSize', 15)
ylim([-4, 4])
xlim([-4, 4])

title('Nyquist diagram (zoomed)')

ax = subplot(2,2,[3,4]), cla
margin(L_cz)
grid on, hold on
% bode(1/(1+L_cz))
% [mx, wmx] = norm(1/(1+L_cz), Inf)
% plot(ax, wmx, 20*log10(mx), 'x')

axs = get(gcf, 'Children');
ylim(axs(3), [-50, 35]);
ylim(axs(2), [-360, 720]);

% subplot(2,2,4)
% rlocus(L_cz)


%%
t = [0:0.01:2*pi]';
x = cos(t);
y = sin(t);


g = zpk([-1, 0.25], [1+j, 1-j], 5);

k = 1;

[re, im] = nyquist(k*g);
re = re(:);
im = im(:);
figure(4); clf
subplot(2,2,2);
nyquist(k*g)

subplot(2,2,1)
plot(re, im)
hold on
plot(re, -im, '--')
plot(x,y, ':k')
grid on

subplot(2,2,3)
rlocus(k*g)

subplot(2,2,4)
bode(k*g)
grid  on
%%

t = [0:0.01:2*pi]';
x = cos(t);
y = sin(t);


g = zpk([0.25], [-.25+j, -.25-j], .75);
% g = zpk([-.025+1.2j, -.025-1.2j], [-.1+j, -.1-j, 0 ], 2);
k = 1;

[re, im] = nyquist(k*g);
re = re(:);
im = im(:);
figure(4); clf

subplot(1,2,1)
margin(g)

subplot(1,2,2);
nyquist(g)

isstable(feedback(g, 1))

figure(5)
rlocus(g)
%%
plot(re, im)
hold on
plot(re, -im, '--')
plot(x,y, ':k')
grid on

subplot(2,2,3)
rlocus(k*g)

subplot(2,2,4)
bode(k*g)
grid  on

margin(g)
grid on

isstable(feedback(g,1))