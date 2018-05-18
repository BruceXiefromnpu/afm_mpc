classdef CanonCntrlParams_01
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
    
    function self = CanonCntrlParams_01(sys)
    % if we make the last one 0.4, the high freq poles with
    % zeta ~=0.5 actually lose damping.
    % self.zeta_s = [.8, .7, .4, .4 .4];
      self.zeta_s = [.8, .7, .4, .4 .542];
      self.gam_s = [1., 1., 1., 1., 1];
      self.rad = 0.25;
      self.pint = 0.8;
      
      [wp_real_x, wz_real_x] = w_zp_real(sys);
      rho_1 = wz_real_x(1)/wp_real_x(1);
      % self.rho_s = [rho_1, 1., 1];
      self.rho_s = [1.02, 1., 1];
      
    end
    
    
  end
  
        
        
end

    
    