% charDes = getCharDes(sys, gam_s, p_int, zeta_des_s, rho_s, rad);
% This is a simple algo to get desired pole locations for a discrete SS
% representation. The algo leaves all real poles and poles at z=0 untouched
% and moves the rest so that the have damping of zeta = 0.7, and increases 
% wn by a factor gam
% 
% sys is a discrete state space system.
% 
% gam_wns is a scalar specifying amount to increase nat frequency by:
% wn <-- wn*gam
%
% rad is not needed if the system has no delay (poles at z=0).


function charDes = getCharDes_const_sig(sys, p_int, cmplx_rad, rho_s, rad);

[wn, z, pls] = damp(sys);
Ts = sys.Ts;
kcp = find(imag(pls) > 0);
kcm  = find(imag(pls) < 0);
% kc_r  = find(imag(pls) == 0);


% if ~exist('zds', 'var')
%     zd = .3*ones(length(kcp));
% end

if  ~exist('rad', 'var')
    rad = 0;
end


k_real = find(imag(pls)==0);
Nd     = length(find(pls(k_real) ==0) );

delPls = getDelayPls(Nd, rad);
for i = 1:length(kcp)
    p_i = pls(kcp(i));
    theta = atan2(imag(p_i), real(p_i));
    p_i_place = cmplx_rad*cos(theta) + cmplx_rad*sin(theta)*1j;
%     wd = wn(kcp(i))*sqrt(1-zd^2)*gam_s(i);
%     sig = zd*wn(kcp(i))*gam_s(i);
%     s = -sig + j*wd;
    charDes(kcp(i)) = p_i_place;
    charDes(kcm(i)) = conj(p_i_place);
end
% keyboard
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


function delPls = getDelayPls(Nd, r)
% given radius r and Nd, generate evenly spaced circle of delay poles, with one
% always left at the origin.

  Nd = Nd - 1; % always leave one at the origin.
  % Put them at the Nd-1 roots of unity.

  k_s = [0:1:Nd-1]; % Take every other one.

  % We have to ensure complex-conjugacy, so we round to the 10th decimal. This
  % is fairly heurstic, but the error from conjugacy seems to on the order of
  % 10^-15 or so.
  delPls = round(r*exp(2j*pi*k_s/Nd), 10);
  delPls = [delPls.'; 0];
  

% ------ Old Method -----. This will only leave a pole at z=0 for Nd odd.

%     Nd_even = floor(Nd/2)+1;
%     div = pi/Nd_even;
%
%     thets = div:div:pi-div;
%
%     rs    = 0*thets + r;
%     p1 = rs.*cos(thets) + rs.*sin(thets)*1i;
%
%     delPls = [p1, conj(p1)];
%
%     if floor(Nd/2) ~= Nd/2
%         delPls = [delPls, 0];
%     end
end

end