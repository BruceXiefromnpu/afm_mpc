clear;
clc;

addpath('../afm_mpc_journal/functions')
addpath('../afm_mpc_journal/functions/canon')
% ---------- Load Parametric Models  -----------
[plants, frf_data] = CanonPlants.plants_ns14();
sys = ss(plants.Gvib);

Ts = sys.Ts;
sys.InputDelay = 9;
sys = absorbDelay(sys);
wnyq = (pi/Ts);
freqs = frf_data.freqs_Hz;


% pint = 0.9;
% Px = getCharDes(sys_recyc, [1, 1], pint, [0.7], [2], 0.25);
can_cntrl = CanonCntrlParams_ns14();
[Q, R0, S1] = build_control(sys, can_cntrl);
can_cntrl.rad = 0.1;
can_cntrl.rho_s = [1, 1];
can_cntrl.gam_s = [1. 1. 1.5 1 1 1];
can_cntrl.zeta_s = [0.5, 0.5, 0.6, 0.4 0.4, 0.4]

[Q, R0, S] = build_control(sys, can_cntrl);
S-S1
% [chat, dhat] = place_zeros(sys_recyc, Px);
% Q = chat'*chat;
% S = chat'*dhat;

%
F1 = figure(1); clf
F2 = figure(2); clf

gam_s = [0.0001, 1, 3];
H_Sens = gobjects(1, length(gam_s));
H_hyd = gobjects(1, length(gam_s));
H_distresp = gobjects(1, length(gam_s));
H_stepresp = gobjects(1, length(gam_s));
CS = {'k', 'g', 'r', 'k', 'm'};

frfBode(sys, freqs, F1, 'Hz', ':k');

for k = 1:length(gam_s)
  gam = gam_s(k);
  R = R0 + gam;

  K_lqr = dlqr(sys.a, sys.b, Q, R, S);
  Kx = K_lqr;

  
  Nbar = SSTools.getNbar(sys, K_lqr);
  Nx = SSTools.getNxNu(sys);

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

  AA_ = [A_tilde, -(sys.b*Nbar+Lx*Cd);
         -Ld*sys.c,    1-Ld*Cd];

  BB_ = [sys.b*Nbar; 0];
  LL_ = [Lx; Ld];

  K_c = [K_lqr, Nbar];
  D1 = ss(AA_, LL_, K_c, 0, Ts);
  D2 = ss(AA_, BB_, -K_c, Nbar, Ts);

  % % Estimated state feedback loop gain:
  loop = sys*D1;
  Sens = sys/(1+loop);

  H_yr = ((sys*D2*Sens));
  H_yd = minreal( sys*Sens);

  %H_hyd(k) = frfBode(H_yd, omegas, F1, ['-', CS{k}], 'rad');
  
  H_Sens(k) = frfBode(Sens, freqs, F1, 'Hz', ['-', CS{k}]);
  hold on
  %H_hyd(k).DisplayName = sprintf('$H_{yd}=G\\mathcal{S}$, $\\gamma=%.4f$', gam);
  H_Sens(k).DisplayName = sprintf('$\\mathcal{S}$, $\\gamma=%.4f$', gam);
                                  

  % leg = legend([H_hyd(1:k), H_Sens(1:k)]);
  leg = legend([H_Sens(1:k)]);
  set(leg, 'FontSize', 16, 'Location', 'SouthEast')

  PLANT = sys;

  % figure(9)
  % pzplot(S, H_yd)
  % title('H_yd, Sensitivity')

  ref = 0;
  dist = 1;
  x0 = PLANT.b*0;
  x0_obs = sys_obs.b*0;
  trun = 400*Ts;
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
  H_distresp(k) = plot(y_dist_sim.Time, y_dist_sim.Data, CS{k});
  H_distresp(k).DisplayName = sprintf('$\\gamma=%.4g$', gam);
  hold on
  %plot(y_dist_sim.Time, y_dist_TF, '--k')
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
