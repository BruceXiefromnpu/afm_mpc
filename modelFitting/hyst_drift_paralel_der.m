

function [err, Derr, u_drift] = hyst_drift_paralel_der(theta, u, y, t, Nhyst, Nsat, nsd, Gvib, r, d)
  
  w_hyst = theta(1:Nhyst);
  w_sat = theta(Nhyst+1:Nhyst+Nsat);
  
  theta_gd = theta(Nhyst+1+Nsat:end);
  %theta_gd = theta(N_hyst+1:end);
  
  lams = (theta_gd(1:nsd));
  c = theta_gd(nsd+1:end-1)';
  b = c'*0+1;
  Dfeed = theta_gd(end);

  if nargout >1
    [u_hyst, DSH_dw_ws] = PIHyst.hyst_play_sat_op(u, r, w_hyst, d, w_sat, w_hyst*0);
    [u_drift, dGd_dlam_c] = recurse_drift(u, lams, c, b, Dfeed);
    Derr = [DSH_dw_ws, dGd_dlam_c];
    %ZeroStateResp = CondenCon.zero_state_output_resp(ss(Gvib), length(u)-1, 0);
    u_effec = (u_drift + u_hyst);
    
    % Mutliplication by the zero state resp matrix is equivivalent to running
    % each of u_effec and each column of the Jacobian through lsim, which is
    % more memory efficient.
    
    %y_fit = ZeroStateResp*u_effec;
    %Derr = ZeroStateResp*Derr;
    %y_fit = u_effec;

    y_fit = lsim(Gvib, u_effec, t);
    for k=1:size(Derr, 2) % across columns
      Derr(:,k) = lsim(Gvib, Derr(:,k), t);
    end
    
  else
    [u_hyst] = PIHyst.hyst_play_sat_op(u, r, w_hyst, d, w_sat, w_hyst*0);
    gd = ss(diag(lams), b, c, Dfeed, Gvib.Ts);
    u_drift = lsim(gd, u, t);
    u_effec = u_drift + u_hyst;
    y_fit = lsim(Gvib, u_effec, t);
    %y_fit = u_effec;
  end
  
  err = y_fit - y;
end

function [y_drift, dGd_dlam_cd] = recurse_drift(u, lams, c, b, d)
  
  %x_mat = zeros(length(lams), length(u));
  y_drift = zeros(length(u), 1);
  xk = lams(:)*0;
  xk_min1 = xk;
  dxk_dlam = lams(:)*0;
  dGd_dlam_cd = zeros(length(u), length(lams)*2+1);
  int_vec = (length(u):-1:1);
  lam_pow = lams.^(int_vec-1);
  
  for k=1:length(u)
    
    y_drift(k) = c*xk + d*u(k);
    
    dxk_dlam = xk_min1 + lams.*dxk_dlam;
    %int_vec2 = (k-1:-1:1);
%     if k>1
%       %dxk_dlam2 = int_vec2.*(lams.^(int_vec2-1))*u(1:k-1);
%       dxk_dlam = int_vec(end-k+2:end).*lam_pow(:,end-k+2:end)*u(1:k-1);
%      %keyboard
%     else
%       dxk_dlam = 0*lams(:);
%     end
    dGd_dlam_cd(k, :) = [(dxk_dlam)'.*c , xk', u(k)];
    
    xk_min1 = xk;
    xk = xk.*lams(:) + b*u(k);
  end
  
end