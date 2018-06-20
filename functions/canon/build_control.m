function [Q, R0, S, P_x] = build_control(sys_recyc, can_cntrl)
% [Q, R0, S] = build_control(sys_recyc, can_cntrl)
% Construct the quadratic cost matrices using the root-locus
% technique (fictitious zeros). 
% 
% Inputs
% ------
%  sys_recyc : the deltaUk ss system
%  can_control : an instance of CanonControlParams_XX.m
  
  
  if sys_recyc.InputDelay ~= 0 || sys_recyc.IOdelay ~=0
    error(['All Delays of sys_obs must be part of the states, not ' ...
           'included as a property']);
  end
  
  
  P_x  = getCharDes(sys_recyc, can_cntrl.gam_s, can_cntrl.pint,...
                               can_cntrl.zeta_s, can_cntrl.rho_s, can_cntrl.rad);
  [Chat, Dhat] = place_zeros(sys_recyc, P_x);
  
  Q = Chat'*Chat;
  S = Chat'*Dhat;
  R0 = Dhat'*Dhat;
  

end
