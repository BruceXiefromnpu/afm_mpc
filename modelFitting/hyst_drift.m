

function [err, Derr] = hyst_drift(theta, u, y, t, N_hyst, Nsat, r, d)
  
  w_hyst = theta(1:N_hyst);
  w_sat = theta(N_hyst+1:N_hyst+Nsat);
  
  [u_hyst, Derr] = PIHyst.hyst_play_sat_op(u, r, w_hyst, d, w_sat, w_hyst*0);

  
%   try
%   y_sim = lsim(Gvib, u_effective, t);
%   catch
%     keyboard
%   end
  err = u_hyst - y;
end