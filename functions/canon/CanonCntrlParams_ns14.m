classdef CanonCntrlParams_ns14
% Constructs the standard pole-placement data for the AFM.
% zeta_s = [.8, .7, .4, .4 .542];
% gams = [1., 1., 1., 1., 1];
% rad = 0.25
% pint = 0.8
% rhos = [1.02, 1,1]
  
  properties
    rho_s;
    gam_s;
    zeta_s;
    rad;
    pint;
    
  end

  methods
    
    function self = CanonCntrlParams_ns14(sys)
      
    % if we make the last one 0.4, the high freq poles with
    % zeta ~=0.5 actually lose damping.
      self.zeta_s = [.85, 0.85, .7, .4, .4 .4];
      self.gam_s = [1., 1., 1., 1., 1., 1.]; % puts it at the extant zero
      self.rad = 0.25;
      self.pint = 0.8;
      
      %[wp_real_x, wz_real_x] = w_zp_real(sys);
%       rho_1 = wz_real_x(1)/wp_real_x(1);
      % self.rho_s = [rho_1, 1., 1];
      self.rho_s = [1.2, 1., 1];
      
    end
    
    function self = aggressive_params(self)
      self.zeta_s = [.7, 0.7, .6, .4, .4 .4];
      self.gam_s = [1.25, 1.2, 1.3, 1., 1., 1.]; % puts it at the extant zero
      self.rad = 0.25;
      self.pint = 0.8;
      self.rho_s = [2.0, 1., 1];
    end
    
  end
  
        
        
end

    
    