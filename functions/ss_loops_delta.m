
function [Sens, Hyd, Hyr, Hyeta, Loop] =ss_loops_delta(sys, sys_recyc, KxKu, Lx)
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


  A_tilde = sys.a - sys.b*Kx - Lx*sys.c;

  AA_ = [A_tilde, sys.b-sys.b*Ku;
       -Kx,     1-Ku,            ];


  BB_ = [sys.b*Nbar; Nbar];
  LL_ = [Lx; 0];

  K_c = KxKu;
  K_c(end) = K_c(end)-1;
  D1 = ss(AA_, LL_, K_c, 0, Ts);
  D2 = ss(AA_, BB_, -K_c, Nbar, Ts);

  % % Estimated state feedback loop gain:
  Loop = sys*D1;
  Sens = 1/(1+Loop);

  Hyr = ((sys*D2*Sens));
  Hyd = minreal( sys*Sens);
  
  Hyeta = -minreal(feedback(sys*D1, 1));

end