
function [Sens, Hyd, Hyr, Hyeta, Loop] =ss_loops(sys, Kx, Lx)
% [Sens, Hyd, Hyr] =ss_loops_delta_dist(sys,  Kx, Lx)
% 
% For the closed loop vanilla s.s. system, constructs the closed loop transfer
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
  
  Nbar = SSTools.getNbar(sys, Kx);

  A_tilde = sys.a - sys.b*Kx - Lx*sys.c;

  D1 = ss(A_tilde, Lx, Kx, 0, Ts);
  D2 = ss(A_tilde, sys.b, -Kx, 1, Ts);

  % % Estimated state feedback loop gain:
  Loop = minreal(zpk(sys)*zpk(D1));
  Sens = minreal(1/(1+Loop));

  Hyr = Nbar*zpk(sys)*zpk(D2)*zpk(Sens);
  Hyr = ss(minreal(Hyr));
  
  Hyd = minreal(sys*Sens);
  Hyeta = -minreal(feedback(sys*D1, 1));

  
end