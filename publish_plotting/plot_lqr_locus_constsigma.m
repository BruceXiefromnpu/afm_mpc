% This script plots the LQR based root locus for the "lowgain"
% version of things.

clear, clc
saveon = true;
plants = CanonPlants.plants_ns14();
G = plants.SYS;
G_recyc = plants.sys_recyc;
Ts = G_recyc.Ts;

const_sig = false;
if const_sig
  figbase = 30;
  pmgm_figfile = 'PMGM_vs_gamma_constsig_0p9.svg';
  sens_ts_figfile = 'GainS_TS_vs_gamma_constsig_0p9.svg';
  lqr_locus_figfile = 'lqr_locus_constsig_0p9.svg';
  cmplx_rad = 0.9;
   [Q1, R0, S1, P_x] = build_control_constsigma(G_recyc, cmplx_rad);
else
  figbase = 0;
  pmgm_figfile = 'PMGM_vs_gamma_setzeta.svg';
  sens_ts_figfile = 'GainS_TS_vs_gamma_setzeta.svg';
  lqr_locus_figfile = 'lqr_locus_setzeta.svg';
  
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




f3 = figure(3+figbase); clf
if const_sig
  t = (0:.01:pi);
  x = cmplx_rad*sin(t);
  y = cmplx_rad*cos(t);
  plot(real(P_x), imag(P_x), 'ob')
  hold on
  plot(x, y, 'k')
end

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
  saveas(f3, fullfile(PATHS.jfig, lqr_locus_figfile))
end

Qw = plants.SYS.b*plants.SYS.b'*50;
Lx = G.a*dlqr(plants.SYS.a', plants.SYS.c', Qw, 1)';
p_int_d = 0.7;
[LxLd, G_obsDist, Ident_obs, C_ydist] = DistEst.output_dist_est(G, Lx, p_int_d);

% Estimator


gams = logspace(log10(0.1), log10(300), 100);

[GM_s, PM_s, Sens_gain, TS_s] = gmpm_vs_gam_recyc(G, G_recyc, G_obsDist, Q1, R0, S1, LxLd, gams);


% f4 = figure(4+figbase); clf
clc
width = 3.45;
height = 3;
f4 = mkfig(4+figbase, width, height);
clf
ax1 = axes('Position', [0.1100 0.1300 0.7750 0.8150], 'YAxisLocation', 'left');

plot_gmpm(gams, GM_s, PM_s, ax1, f4);



% f5 = figure(5+figbase); clf
f5 = mkfig(5+figbase, width, height);
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
   saveas(f4, fullfile(PATHS.jfig, pmgm_figfile))
   saveas(f5, fullfile(PATHS.jfig, sens_ts_figfile))
end


%%
  K_lqr = dlqr(G_recyc.a, G_recyc.b, Q1, R0+107, S1);
  [Sens, Hyd, Hyr, Hyeta, Loop] = ss_loops_delta_dist(plants.SYS, G_recyc, G_obsDist, K_lqr, LxLd);
  figure, step(Hyr, Hyd)
  








