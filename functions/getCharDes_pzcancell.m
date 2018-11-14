function charDes = getCharDes_pzcancell(sys, cancel_cmplx_idx, gam_s, p_int, zet_des_s, rho_s, rad)
  
  p = pole(sys);
  z = tzero(sys);
  
  p = sort_by_w(p);
  z = sort_by_w(z);
  
  
  z_real = z(imag(z) == 0);
  p_real = z(imag(z) == 0);
  
  z_cmplx = z(imag(z) ~= 0);
  p_cmplx = p(imag(p) ~= 0);
  
  
  idx_cur = 1;
  
  Ts = sys.Ts;
  kcp = find(imag(p) > 0);
  kcm  = find(imag(p) < 0);

  wns = abs(log(p_cmplx)/Ts);
  sigs = -real(log(p_cmplx)/Ts);
  zets = sigs./wns;
  
  
  
  
  for i = 1:length(cancel_cmplx_idx)
    z_idx = kcp(cancel_cmplx_idx(i));
    charDes(kcp(i)) = z_cmplx(z_idx);
    charDes(kcm(i)) = conj(z_cmplx(z_idx));
        
  end

  wns_new = wns.*gam_s(1:length(wns));
  zets_new = zet_des_s(1:length(zets));
  pdes_cmplx = wns_new.*zets_new +
  
%   delPls = getDelayPls(Nd, rad);
%   for i = 1:length(kcp)
%       zd = zeta_des_s(i);
% 
%       wd = wn(kcp(i))*sqrt(1-zd^2)*gam_s(i);
%       sig = zd*wn(kcp(i))*gam_s(i);
%       s = -sig + j*wd;
% 
%       charDes(kcp(i)) = exp(Ts*s);
%       charDes(kcm(i)) = exp(Ts*conj(s));
% 
%   end

  jj = 1;
  delPls_iter = 1;

  for i=1:length(k_real)
    if pls(k_real(i)) ~= 1 && pls(k_real(i))~=0 || Ts==0
      % Move frequency of real polse by rho_s.
      if Ts ~=0
        s = -wn(k_real(i))*rho_s(jj);
        charDes(k_real(i)) = exp(Ts*s);
        jj = jj+1;
      elseif Ts ==0
        charDes(k_real(i)) = pls(k_real(i))*rho_s(jj);
        jj = jj+1;
      end

    elseif  pls(k_real(i)) == 1 & Ts ~= 0 % place augmented integrator pole
      charDes(k_real(i)) = p_int;

    else % place delay poles.
      %                  charDes(k_real(i)) = pls(k_real(i));
      charDes(k_real(i)) = delPls(delPls_iter);
      delPls_iter = delPls_iter+1;
    end



  end
  
end

