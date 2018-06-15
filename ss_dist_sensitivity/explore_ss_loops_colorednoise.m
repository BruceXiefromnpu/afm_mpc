clear;
clc;



Ts = 0.01;
wz1 = 5
wp1 = 5.1
g2 = tf([1 2*0.01*wz1, wz1^2], [1, 2*0.01*wp1, wp1^2]);
g3 = tf(8, [1 8]);
sys = c2d(ss(tf([144], [1, 2*12*0.01, 144])), Ts);
sys.InputDelay =0;
sys = absorbDelay(sys);
% figure(5)
% bode(sys)

wnyq = (pi/Ts);
omegas1 = logspace(log10(0.005), log10(wnyq), 200);
[~,~, omegas] = bode(sys);
omegas = unique([omegas1(:); omegas(:)] );


pint = 0.8;
Px = getCharDes(sys, [1, 1], pint, [0.85, 0.85], [2], .25);
[chat, dhat] = place_zeros(sys, Px);
Q = chat'*chat;
S = chat'*dhat;
R0 = dhat^2;

gam_s = [.0001, 3, 15]


addpath('../afm_mpc_journal/functions')
addpath('../afm_mpc_journal/functions/canon')
% ---------- Load Parametric Models  -----------

clc
Ts = sys.Ts;
sys.InputDelay = 0;
sys = absorbDelay(sys);
wnyq = (pi/Ts);
freqs = omegas/2/pi;

F1 = figure(1); clf
F2 = figure(2); clf

H_Sens = gobjects(1, length(gam_s));
H_hyd = gobjects(1, length(gam_s));
H_distresp = gobjects(1, length(gam_s));
H_stepresp = gobjects(1, length(gam_s));

plotstyle = {'color', 'k', 'LineStyle', ':'};
frfBode(sys, freqs, F1, plotstyle, 'Hz');

mp = colormap(gca, 'jet');
x = 1:size(mp,1);
xq = linspace(1, (size(mp,1)), length(gam_s));
mp_fine = zeros(length(xq), 3);
mp_fine(:,1) = interp1(x, mp(:,1), xq);
mp_fine(:,2) = interp1(x, mp(:,2), xq);
mp_fine(:,3) = interp1(x, mp(:,3), xq);




Gu = tf([1, 2*pi*.01], [1, 2*pi*0.1]);
Gu = ss(c2d(Gu, Ts))*0.01;
frfBode(Gu, freqs, F1, plotstyle, 'Hz');
qw = eye(length(Gu.a));
abar = [sys.a, sys.b*Gu.c; 
        sys.c*0, Gu.a];
bbar = [sys.b; 0];
cbar = [sys.c, 0];
ebar = [sys.b*Gu.d; Gu.b];

sys_bar = ss(abar, bbar, cbar, 0, Ts);

Qwbar = ebar*qw*ebar'*5000;

% Qw = sys.b*sys.b'*50;
Lx = dlqr(sys_bar.a', sys_bar.c', Qwbar, 1)';
% sys_obs = sys_bar;
% sys_bar = sys;
% Lx = dlqr(sys_bar.a', sys_bar.c', sys_bar.b*sys_bar.b'*100, 1)';
% L = Lx;


for k = 1:2 %length(gam_s)
  gam = gam_s(k);
  R = R0 + gam;

  Kx = dlqr(sys.a, sys.b, Q, R, S);

  Nbar = SSTools.getNbar(sys, Kx);
  Nx = SSTools.getNxNu(sys);

  p_int_d = 0.8;
  Lx = dlqr(sys_bar.a', sys_bar.c', Qwbar, 1)';
  [L, sys_obs, Ident_obs, C_ydist] = DistEst.output_dist_est(sys_bar, Lx, p_int_d);

  Ld = L(end);
  Lx = L(1:end-1);
  Cd = sys_obs.c(end);

  % Create ss systems.
  A_tilde = sys_bar.a - sys_bar.b*[Kx, 0] - Lx*sys_bar.c;

  AA_ = [A_tilde,  -(sys_bar.b*Nbar+Lx*Cd);
         -Ld*sys_bar.c,           1-Ld*Cd];

  BB_ = [sys_bar.b*Nbar; 0];
  LL_ = [Lx; Ld];
  K_c = [Kx, 0, Nbar];
  D1 = ss(AA_, LL_, K_c, 0, Ts);
  D2 = ss(AA_, BB_, -K_c, Nbar, Ts);

%   % % Estimated state feedback loop gain:
  loop = sys*D1;
  Sens = 1/(1+loop);
 
  H_yr = ((sys*D2*Sens));
  H_yd = minreal( sys*Sens);
 
%   %H_hyd(k) = frfBode(H_yd, omegas, F1, ['-', CS{k}], 'Hz');
  clr = mp_fine(k,:);
  plot_style = {'color', clr};
  H_Sens(k) = frfBode(Sens, freqs, F1, plot_style, 'Hz');
%   hold on
%   %H_hyd(k).DisplayName = sprintf('$H_{yd}=G\\mathcal{S}$, $\\gamma=%.4f$', gam);
%   H_Sens(k).DisplayName = sprintf('$\\mathcal{S}$, $\\gamma=%.4f$', gam);
% 
  leg = legend([H_Sens(1:k)]);
  set(leg, 'FontSize', 16, 'Location', 'SouthEast')
% 
  PLANT = sys;

  ref = 0;
  dist = 1;
  x0 = PLANT.b*0;
  x0_obs = sys_obs.b*0;
  trun = 300*Ts;
  K_lqr = [Kx, 0];
  sim('AFMss_dist_obshas_uk')
  y_dist_sim = Y;
  u = y_dist_sim.Time*0+ref;
  d = 0*u + dist;
  y_dist_TF = lsim(H_yd, d, y_dist_sim.Time);

  ref = 1;
  dist = 0;
  sim('AFMss_dist_obshas_uk')
  y_ref_sim = Y;
  u = y_dist_sim.Time*0+ref;
  y_ref_TF = lsim(H_yr, u, y_dist_sim.Time);


  figure(F2);
  subplot(2,1,1)
  H_distresp(k) = plot(y_dist_sim.Time, y_dist_sim.Data, plot_style{:});
  H_distresp(k).DisplayName = sprintf('$\\gamma=%.4g$', gam);
  hold on
%   plot(y_dist_sim.Time, y_dist_TF, '--')
  title('Disturbance Response')
  grid on
  
  subplot(2,1,2)
  H_stepresp(k) = plot(y_ref_sim.Time, y_ref_sim.Data, plot_style{:});
  H_distresp(k).DisplayName = sprintf('$\\gamma=%.4g$', gam);
  hold on
  %plot(y_ref_sim.Time, y_ref_TF, '--')
  title('Step Response')
  legend(H_distresp(1:k))
  grid on
end
return





%%
gam_s = logspace(log10(0.0001), log10(25),25);
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
  
  
dcgs = gam_s*0;
dcg_num = gam_s*0;
dcg_den = gam_s*0;
I = zpk([], [1], 1, Ts);

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
  
  g_tmp = minreal(I*zpk(Sens));
  [zz, pp, kk] = zpkdata(g_tmp);
  gnum = zpk(zz, [], 1, Ts);
  gden = zpk([], pp, 1, Ts);
  dcg_num(k) = evalfr(gnum, 1);
  dcg_den(k) = evalfr(gden, 1);
  
  dcgs(k) = dcgain(minreal(I*zpk(Sens)));
  if k == 1
    s1 = Sens;
  elseif k == length(gam_s)
    s2 = Sens;
  end

  H_yr = ((zpk(sys)*zpk(D2)*zpk(Sens)));
  H_yd = minreal( zpk(sys)*zpk(Sens));

  %H_hyd(k) = frfBode(H_yd, omegas, F1, ['-', CS{k}], 'Hz');
  plot_style = {'color', clr};
  H_Sens(k) = frfBode(Sens, freqs, F3, plot_style, 'Hz');
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

%%
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
%%
% ------------------------------------------------------------------------
CS = {'k', 'g', 'r', 'k', 'm'};
for k = 1:length(gam_s)
  gam = gam_s(k);
  R = dhat^2+gam;

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
  Qw = sys.b*sys.b'*10;
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

%   H_hyd(k) = frfBode(H_yd, omegas, F1, ['-', CS{k}], 'rad');
  
  H_Sens(k) = frfBode(Sens, omegas, F1, ['--', CS{k}], 'rad');
hold on
%   H_hyd(k).DisplayName = sprintf('$H_{yd}=G\\mathcal{S}$, $\\gamma=%.4f$', gam);
  H_Sens(k).DisplayName = sprintf('$\\mathcal{S}$, $\\gamma=%.4f$', gam);
                                  

%   leg = legend([H_hyd(1:k), H_Sens(1:k)]);
  leg = legend(H_Sens(1:k));
  set(leg, 'FontSize', 16, 'Location', 'SouthEast')

  PLANT = sys;

  % figure(9)
  % pzplot(S, H_yd)
  % title('H_yd, Sensitivity')

  ref = 0;
  dist = 1;
  x0 = PLANT.b*0;
  x0_obs = sys_obs.b*0;
  trun = 6;
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
  subplot(2,1,1)
  H_distresp(k) = plot(y_dist_sim.Time, y_dist_sim.Data, CS{k});
  H_distresp(k).DisplayName = sprintf('$\\gamma=%.4g$', gam);
  hold on
  %plot(y_dist_sim.Time, y_dist_TF, '--')
  title('Disturbance Response')

  subplot(2,1,2)
  H_stepresp(k) = plot(y_ref_sim.Time, y_ref_sim.Data, CS{k});
  H_distresp(k).DisplayName = sprintf('$\\gamma=%.4g$', gam);
  hold on
  %plot(y_ref_sim.Time, y_ref_TF, '--')
  title('Step Response')
  legend(H_distresp(1:k))

end
return


pole(H)


loop = minreal(sys*D1);
M = [loop.a, loop.b;
     loop.c, 1];

P = eye(size(M,1));
P(end,end) = 0;
pp = eig(M, P)
