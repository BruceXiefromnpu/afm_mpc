
function [Sens, Hyd, Hyr, Loop] =ss_loops_delta_dist_frf(G_frf, omegas, sys, sys_recyc, sys_obs, KxKu, LxLd)
% [Sens, Hyd, Hyr] =ss_loops_delta_dist(sys, sys_recyc, sys_obs, KxKu, LxLd)
% 
% For the closed loop s.s. system with deltaU augmentation and disturbance
% estimation (input disturbance only), constructs the closed loop transfer
% functions:
%
% Sens = 1/(1 + Loop)
% Loop = G(z)* K_c(zI - A_tilde)^-1 * LL_
%
% Hyr = G(z)[Nbar - K_c(zI - A_tilde)^-1 BB_ ]
%       --------------------------------------
%        1 + G(z) * K_c(zI - A_tilde)^-1 * LL_
%
% Hyd =                 G(z)
%       --------------------------------------
%        1 + G(z) * K_c(zI - A_tilde)^-1 * LL_
%
  Ts = sys.Ts;
  
  Kx = KxKu(1:end-1);
  Ku = KxKu(end);
  Nbar = SSTools.getNbar(sys_recyc, KxKu);

  Ld = LxLd(end);
  Lx = LxLd(1:end-1);
  Cd = sys_obs.c(end);

  A_tilde = sys.a - sys.b*Kx - Lx*sys.c;

  AA_ = [A_tilde, sys.b-sys.b*Ku, -(sys.b*Nbar+Lx*Cd);
       -Kx,     1-Ku,             -Nbar;
       -Ld*sys.c,    0,              1-Ld*Cd];

  BB_ = [sys.b*Nbar; Nbar; 0];
  LL_ = [Lx; 0; Ld];

  K_c = [KxKu, Nbar];
  K_c(end-1) = K_c(end-1)-1;
  D1 = ss(AA_, LL_, K_c, 0, Ts);
  D2 = ss(AA_, BB_, -K_c, Nbar, Ts);

  D1_frf = squeeze(freqresp(D1, omegas));
  D2_frf = squeeze(freqresp(D2, omegas));
  
  % % Estimated state feedback loop gain:
  Loop = G_frf.*D1_frf;
  Sens = 1./(1+Loop);

  Hyr = G_frf.*D2_frf.*Sens;
  Hyd = G_frf.*Sens;

end