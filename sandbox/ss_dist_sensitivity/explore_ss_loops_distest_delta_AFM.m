clear;
clc;

addpath('../functions')
addpath('../functions/canon')
% ---------- Load Parametric Models  -----------
[plants, frf_data] = CanonPlants.plants_ns14();
sys = ss(plants.Gvib);

Ts = sys.Ts;
sys.InputDelay = 9;
sys = absorbDelay(sys);
wnyq = (pi/Ts);
freqs = frf_data.freqs_Hz;

sys_recyc = SSTools.deltaUkSys(sys);
% sys_recyc.b(end-1) = 0;
% pint = 0.9;
% Px = getCharDes(sys_recyc, [1, 1], pint, [0.7], [2], 0.25);
% can_cntrl.gam_s = [1.5 1.5 1.5 1 1 1];
% can_cntrl.zeta_s = [0.5, 0.6, 0.6, 0.6 0.6, 0.6]
% P_x  = getCharDes(sys_recyc, gam_s, pint, zeta_s, rho_s, rad);

can_cntrl = CanonCntrlParams_ns14();
zeta_s = [.7, 0.77, .7, .4, .4 .4];
gam_s = [1., 1., 1.2, 1., 1., 1.]; % puts it at the extant zero
rad = 0.51;
pint = 0.8;
rho_s = [1, 1., 1];
cmplx_rad = 0.92;
P_x  = getCharDes_const_sig(sys_recyc, pint, cmplx_rad, rho_s, rad);
% P_x  = getCharDes(sys_recyc, gam_s, pint, zeta_s, rho_s, rad);
[Chat, Dhat] = place_zeros(sys_recyc, P_x);

Q = Chat'*Chat;
S = Chat'*Dhat;
R0 = Dhat'*Dhat;


figure(20); clf
K0 = dlqr(sys_recyc.a, sys_recyc.b, Q, R0+.01, S);
syscl = sys_recyc;
syscl.a = syscl.a - syscl.b*K0;
pzplot(syscl)

theta = 0:0.01:2*pi;
x = cmplx_rad*cos(theta);
y = cmplx_rad*sin(theta);
hold on
plot(x, y, ':k')
plot(real(P_x), imag(P_x), '.', 'MarkerSize', 8)


% [Q, R0, S] = build_control(sys_recyc, can_cntrl);



F1 = figure(1); clf
subplot(2,1,1);
ax1 = gca();
subplot(2,1,2)
ax2 = gca();

F2 = figure(2); clf
ax3 = subplot(2,1,1);
ax4 = subplot(2,1,2);

gam_s = [0.1, 10, 50, 200];
H_Sens = gobjects(1, length(gam_s));
H_hyr = gobjects(1, length(gam_s));
H_distresp = gobjects(1, length(gam_s));
H_stepresp = gobjects(1, length(gam_s));

plotstyle = {'color', 'k', 'LineStyle', ':'};
frfBodeMag(sys, freqs, ax1, 'Hz', plotstyle);

frfBodeMag(sys, freqs, ax2, 'Hz', plotstyle);

% mp = colormap(gca, 'jet');
mp = colormap(gca, 'copper');
x = 1:size(mp,1);
xq = linspace(1, (size(mp,1)), length(gam_s));
mp_fine = zeros(length(xq), 3);
mp_fine(:,1) = interp1(x, mp(:,1), xq);
mp_fine(:,2) = interp1(x, mp(:,2), xq);
mp_fine(:,3) = interp1(x, mp(:,3), xq);


for k = 1:length(gam_s)
  gam = gam_s(k);
  R = R0 + gam;

  Qp = dare(sys_recyc.a, sys_recyc.b, Q, R, S);
  mpc_prob = condensedMPCprob_OA(sys_recyc, 12, Q, Qp, R, S);
  kappa = mpc_prob.kappa;
  
  K_lqr = dlqr(sys_recyc.a, sys_recyc.b, Q, R, S);
  Kx = K_lqr(1:end-1);
  Nbar = SSTools.getNbar(sys_recyc, K_lqr);
  Nx = SSTools.getNxNu(sys_recyc);

  Ku = K_lqr(end);

  Qw = sys.b*sys.b'*100;
  Lx = dlqr(sys.a', sys.c', Qw, 1)';
  p_int_d = 0.8;
  [L, sys_obs, Ident_obs, C_ydist] = DistEst.output_dist_est(sys, Lx, p_int_d);
  Ld = L(end);
  Lx = L(1:end-1);
  Cd = sys_obs.c(end);

  % Create ss systems.
  A_tilde = sys.a - sys.b*Kx - Lx*sys.c;

  AA_ = [A_tilde, sys.b-sys.b*Ku, -(sys.b*Nbar+Lx*Cd);
         -Kx,     1-Ku,             -Nbar;
         -Ld*sys.c,    0,              1-Ld*Cd];

  BB_ = [sys.b*Nbar; Nbar; 0];
  LL_ = [Lx; 0; Ld];

  K_c = [K_lqr, Nbar];
  K_c(end-1) = K_c(end-1)-1;
  D1 = ss(AA_, LL_, K_c, 0, Ts);
  D2 = ss(AA_, BB_, -K_c, Nbar, Ts);

  % % Estimated state feedback loop gain:
  loop = sys*D1;
  Sens = 1/(1+loop);

  H_yr = ((sys*D2*Sens));
  H_yd = minreal( sys*Sens);

  
  
  clr = mp_fine(k,:);
  plot_style = {'color', clr};
  H_hyr(k) = frfBodeMag(H_yr, freqs, ax1, 'Hz', plot_style);
  H_Sens(k) = frfBodeMag(Sens, freqs, ax2, 'Hz', plot_style);
  hold on
%   H_hyd(k).DisplayName = sprintf('$H_{yd}=G\\mathcal{S}$, $\\gamma=%.4f$', gam);
  H_Sens(k).DisplayName = sprintf('$\\mathcal{S}$, $\\gamma=%.4f$, $\\kappa = %.1f$', gam, kappa);
                                  

  % leg = legend([H_hyd(1:k), H_Sens(1:k)]);
  leg = legend([H_Sens(1:k)]);
  set(leg, 'FontSize', 16, 'Location', 'SouthEast')

  PLANT = sys;

  ref = 0;
  dist = 1;
  x0 = PLANT.b*0;
  x0_obs = sys_obs.b*0;
  trun = 400*Ts;
  sim('AFMss_recyc_dist_obshas_uk')
  y_dist_sim = Y;
  u = y_dist_sim.Time*0+ref;
  d = 0*u + dist;
  y_dist_TF = lsim(H_yd, d, y_dist_sim.Time);

  ref = 1;
  dist = 0;
  sim('AFMss_recyc_dist_obshas_uk')
  y_ref_sim = Y;
  u = y_dist_sim.Time*0+ref;
  y_ref_TF = lsim(H_yr, u, y_dist_sim.Time);


  figure(F2);
  H_distresp(k) = plot(ax3, y_dist_sim.Time, y_dist_sim.Data, plot_style{:});
  H_distresp(k).DisplayName = sprintf('$\\gamma=%.4g$, $\\kappa = %.1f$', gam, kappa);
  hold(ax3, 'on')
  %plot(y_dist_sim.Time, y_dist_TF, '--')
  title(ax3, 'Disturbance Response')
  grid(ax3, 'on')

  figure(F2)
  H_stepresp(k) = plot(ax4, y_ref_sim.Time, y_ref_sim.Data, plot_style{:});
  H_distresp(k).DisplayName = sprintf('$\\gamma=%.4g$, $\\kappa = %.1f$', gam, kappa);
  hold(ax4, 'on')
  grid(ax4, 'on')
  %plot(y_ref_sim.Time, y_ref_TF, '--')
  title(ax4, 'Step Response')
  legend(H_distresp(1:k))
  grid(ax4, 'on')
end
return





%%
gam_s = logspace(log10(0.0001), log10(150),25);
% gam_s = linspace((0.0001), (15),25);
% gam_s = [0.0001, 1, 3];
F3 = figure(30); clf;
F9 = figure(90); clf;
mp = colormap(gca, 'jet');
zgrid
F10 = figure(100); clf;
zgrid

mp = colormap(gca, 'jet');
x = 1:size(mp,1);
xq = linspace(1, (size(mp,1)), length(gam_s));
mp_fine = zeros(length(xq), 3);
mp_fine(:,1) = interp1(x, mp(:,1), xq);
mp_fine(:,2) = interp1(x, mp(:,2), xq);
mp_fine(:,3) = interp1(x, mp(:,3), xq);
  
  

j = 1;
for k = 1:length(gam_s)

  clr = mp_fine(k,:);
  
  gam = gam_s(k);
  R = R0 + gam;

  K_lqr = dlqr(sys_recyc.a, sys_recyc.b, Q, R, S);
  % sys_cl = sys_recyc;
  % sys_cl.a = sys_recyc.a - sys_recyc.b*K_lqr;

  % figure(10)
  % hold on
  % pzplot(sys, sys_cl)

  Kx = K_lqr(1:end-1);
  Nbar = SSTools.getNbar(sys_recyc, K_lqr);
  Nx = SSTools.getNxNu(sys_recyc);

  Ku = K_lqr(end);

  sys_obs = sys;
  Qw = sys.b*sys.b'*100;
  Lx = dlqr(sys.a', sys.c', Qw, 1)';
  p_int_d = 0.8;
  [L, sys_obs, Ident_obs, C_ydist] = DistEst.output_dist_est(sys, Lx, p_int_d);
  Ld = L(end);
  Lx = L(1:end-1);
  Cd = sys_obs.c(end);

  % Create ss systems.
  A_tilde = sys.a - sys.b*Kx - Lx*sys.c;

  AA_ = [A_tilde, sys.b-sys.b*Ku, -(sys.b*Nbar+Lx*Cd);
         -Kx,     1-Ku,             -Nbar;
         -Ld*sys.c,    0,              1-Ld*Cd];

  BB_ = [sys.b*Nbar; Nbar; 0];
  LL_ = [Lx; 0; Ld];

  K_c = [K_lqr, Nbar];
  K_c(end-1) = K_c(end-1)-1;
  D1 = ss(AA_, LL_, K_c, 0, Ts);
  D2 = ss(AA_, BB_, -K_c, Nbar, Ts);

  % % Estimated state feedback loop gain:
  loop = sys*D1;
  Sens = 1/(1+loop);

  H_yr = ((zpk(sys)*zpk(D2)*zpk(Sens)));
  H_yd = minreal( zpk(sys)*zpk(Sens));

  %H_hyd(k) = frfBode(H_yd, omegas, F1, 'Hz', ['-', CS{k}]);
  plot_style = {'color', clr};
  H_Sens(k) = frfBode(Sens, freqs, F3, 'Hz', plot_style);
  hold on
  %H_hyd(k).DisplayName = sprintf('$H_{yd}=G\\mathcal{S}$, $\\gamma=%.4f$', gam);
  H_Sens(k).DisplayName = sprintf('$\\mathcal{S}$, $\\gamma=%.4f$', gam);
                                  

  % leg = legend([H_hyd(1:k), H_Sens(1:k)]);
  leg = legend([H_Sens(1:k)]);
  set(leg, 'FontSize', 16, 'Location', 'SouthEast')

  figure(F9)
  SensG = minreal(zpk(Sens)*zpk(sys));
  pg = pole(SensG);
  zg = tzero(SensG);
  plot(real(pg), imag(pg), 'x', 'MarkerSize', 4, 'color', clr);
  plot(real(zg), imag(zg), 'o', 'MarkerSize', 10, 'color', clr);
  title('Sensitivity')
  hold on
  
  figure(F10)
  H_ = minreal(zpk(H_yr));
  ph = pole(H_);
  zh = tzero(H_);
  plot(real(ph), imag(ph), 'x', 'MarkerSize', 4, 'color', clr);
  plot(real(zh), imag(zh), 'o', 'MarkerSize', 4, 'color', clr);
  title('$H_yr$')
  hold on
  
  

end


  tcks_labs = logspace(log10(gam_s(1)), log10(gam_s(end)), 10);
  tcks = linspace(gam_s(1), gam_s(end), 10);
  tcks =  round(tcks, 2, 'significant');
  tcks_labs =  round(tcks_labs, 2, 'significant');
figure(F9);
  caxis([gam_s(1), gam_s(end)])
  C_hand = colorbar('Ticks', tcks, 'TickLabels', tcks_labs, 'Location', 'northoutside'); 
  xlabel('Re');
  ylabel('Im');
  xlim([-1.2, 1.2])
  ylim([-1.2, 1.2])  
  zgrid
figure(F10);
  caxis([gam_s(1), gam_s(end)])
  C_hand = colorbar('Ticks', tcks, 'TickLabels', tcks_labs, 'Location', 'northoutside'); 
  xlabel('Re');
  ylabel('Im');
  zgrid  
  xlim([-1.2, 1.2])
  ylim([-1.2, 1.2])
%
  theta = 0:0.01:2*pi;
  x = cmplx_rad*cos(theta);
  y = cmplx_rad*sin(theta);
  hold on
  plot(x, y, 'k')
  plot(real(P_x), imag(P_x), '.k', 'MarkerSize', 12)