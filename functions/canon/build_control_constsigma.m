function [Q, R0, S, Px] = build_control_constsigma(G_recyc, cmplx_rad)
% [Q, R0, S] = build_control(sys_recyc, can_cntrl)
% Construct the quadratic cost matrices using the root-locus
% technique (fictitious zeros). 
% 
% Inputs
% ------
%  sys_recyc : the deltaUk ss system
%  can_control : an instance of CanonControlParams_XX.m





  p_int = 0.8;  rho_s = [1., 1]; rad = 0.25;
  Px = getCharDes_const_sig(G_recyc, p_int, cmplx_rad, rho_s, rad).';
  
  z = tzero(G_recyc);
  z = sort_by_w(z(imag(z)~=0));
  Px(2:3) = z(1:2);
%   Px(4:5) = z(3:4);
  
  
%   Px(1) = [];
  [Chat, Dhat] = place_zeros(G_recyc, Px);
  Q = Chat'*Chat;
  
  S = Chat'*Dhat;
  R0 = Dhat'*Dhat;
  
%   if sys_recyc.InputDelay ~= 0 || sys_recyc.IOdelay ~=0
%     error(['All Delays of sys_obs must be part of the states, not ' ...
%            'included as a property']);
%   end
%   
%   
%   P_x  = getCharDes(sys_recyc, can_cntrl.gam_s, can_cntrl.pint,...
%                                can_cntrl.zeta_s, can_cntrl.rho_s, can_cntrl.rad);
%   [Chat, Dhat] = place_zeros(sys_recyc, P_x);
%   
%   Q = Chat'*Chat;
%   S = Chat'*Dhat;
%   R0 = Dhat'*Dhat;
  

end
