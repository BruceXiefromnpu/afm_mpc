

function [err, y_sim] = hyst_drift_paralel(theta, u, y, t, N_hyst, Nsat, nsd, Gvib, r, d)
  
  w_hyst = theta(1:N_hyst);
  w_sat = theta(N_hyst+1:N_hyst+Nsat);
  
  theta_gd = theta(N_hyst+1+Nsat:end);
  %theta_gd = theta(N_hyst+1:end);
  
  a = diag(theta_gd(1:nsd));
  c = theta_gd(nsd+1:end-1)';
  b = c'*0+1;
  
  gd = ss(a, b, c, theta_gd(end), Gvib.Ts);
  
  
  u_drft = lsim(gd, u, t);
  
  u_hyst = PIHyst.hyst_play_sat_op(u, r, w_hyst, d, w_sat, w_hyst*0);
  u_effective = u_drft + u_hyst;
%   u_hyst = PIHyst.hyst_play_op(u, r, w_hyst, w_hyst*0);
%   
%   u_effective = PIHyst.sat_op(u_drft + u_hyst, d, w_sat);
  
  try
  y_sim = lsim(Gvib, u_effective, t);
  catch
    keyboard
  end
  err = y - y_sim;
end