% This script plots the LQR based root locus for the "lowgain"
% version of things.

clear, clc
saveon = true;
plants = CanonPlants.plants_ns14();
G_recyc = plants.sys_recyc;
Ts = G_recyc.Ts;

const_sig = true;
if const_sig
  figbase = 30;
else
  figbase = 0;
end
cmplx_rad = 0.9;
G = plants.SYS;
G_recyc = plants.sys_recyc;

if const_sig
  [Q1, R0, S1, P_x] = build_control_constsigma(G_recyc, cmplx_rad);
else
  can_cntrl = CanonCntrlParams_ns14();
%   can_cntrl.zeta_s(1) = 0.7
%   can_cntrl.zeta_s(end-1) = 0.35
%   can_cntrl.pint = 0.9;
  % [Q1, R0, S1] = build_control(G_recyc, can_cntrl);
  P_x  = getCharDes(G_recyc, can_cntrl.gam_s, can_cntrl.pint,...
                               can_cntrl.zeta_s, can_cntrl.rho_s, can_cntrl.rad);
%   z = tzero(G_recyc);
%   z = sort_by_w(z(imag(z)~=0));
%   P_x(2:3) = z(1:2);
  [Chat, Dhat] = place_zeros(G_recyc, P_x);
  Q1 = Chat'*Chat;
  S1 = Chat'*Dhat;
  R0 = Dhat'*Dhat;
end
t = (0:.01:pi);
x = cmplx_rad*sin(t);
y = cmplx_rad*cos(t);

f3 = figure(3+figbase); clf

plot(real(P_x), imag(P_x), 'ob')
hold on
plot(x, y, 'k')

[ax, C_hand] = lqr_locus(G_recyc, Q1, 1, S1, .001, 10000);
C_hand.Label.Interpreter = 'latex';
C_hand.Label.FontSize = 14;
C_hand.Location = 'eastoutside';
C_hand.Position = [0.8406    0.1121    0.0473    0.8715];
ax.Position =  [0.0862 0.1123 0.7369 0.8713];
C_hand.Label.String = '$\gamma$';
C_hand.Label.Units = 'normalized';
C_hand.Label.Rotation = 0;
C_hand.Label.Position = [2.5214    0.5549         0];

xlim([-0.3, 1])
ylim([-0.35, 0.35])
% C_hand.Label.String = '$R_o + \gamma$';

xlabel('Re')
ylabel('Im')


if saveon
  saveas(f3, fullfile(PATHS.jfig, 'lqr_locus_constsig_0p9.svg'))
end

Qw = plants.SYS.b*plants.SYS.b'*50;
Lx = G.a*dlqr(plants.SYS.a', plants.SYS.c', Qw, 1)';
p_int_d = 0.7;
[LxLd, G_obsDist, Ident_obs, C_ydist] = DistEst.output_dist_est(plants.SYS, Lx, p_int_d);



% Estimator


gams = logspace(log10(0.1), log10(1000), 100);

GM_s = [];
PM_s = [];
Sens_gain = [];
Sens_gain2 = [];
TS_s = [];
w_eval = 10;
for k=1:length(gams)
  %
  K_lqr = dlqr(G_recyc.a, G_recyc.b, Q1, R0+gams(k), S1);
  [Sens, Hyd, Hyr, Hyeta, Loop] = ss_loops_delta_dist(plants.SYS, G_recyc, G_obsDist, K_lqr, LxLd);
  
  t = (0:400)'*Ts;
  [y, t] = step(Hyr, t);
  TS_s(k) = settle_time(t, y, 1, 0.01, [], 20);
  
  z = sort(zero(Sens));
  zr = z(imag(z) == 0);
  Int_prox = zpk([], [zr(end)], 1, Ts);
  Sens_gain2(k) = dcgain(minreal(Sens*Int_prox));
  Der = zpk([1], [], 1, Ts);
  der_gain = abs(freqresp(Der, w_eval));
  Sens_gain(k) = abs(freqresp(Sens, w_eval))/der_gain;
  
%   figure(200); hold on; grid on
%   bode(Sens, Der)
%   keyboard
  
  [gm, pm] = margin(Loop);
  
  GM_s(k) = 20*log10(gm);
  PM_s(k) = pm;
end

%%
f4 = figure(4+figbase); clf
yyaxis left

semilogx(gams, GM_s)
hold on
grid on
xlabel('$\gamma$')
ylabel('Gain Margin [dB]')

yyaxis right
plot(gams, PM_s)
ylabel('Phase Margin [deg]')

f5 = figure(5+figbase); clf
yyaxis left
plot(gams, Sens_gain)
hold on
% plot(gams, Sens_gain2, '--')

grid on
xlabel('$\gamma$')
ylabel('Gain of $S(z)$')

yyaxis right
plot(gams, TS_s*1000)
ylabel('Settle-time [ms]')


if saveon
   saveas(f4, fullfile(PATHS.jfig, 'PMGM_vs_gamma_constsig_0p9.svg'))
   saveas(f5, fullfile(PATHS.jfig, 'GainS_TS_vs_gamma_constsig_0p9.svg'))
end


%%
  K_lqr = dlqr(G_recyc.a, G_recyc.b, Q1, R0+107, S1);
  [Sens, Hyd, Hyr, Hyeta, Loop] = ss_loops_delta_dist(plants.SYS, G_recyc, G_obsDist, K_lqr, LxLd);
  figure, step(Hyr, Hyd)
  








