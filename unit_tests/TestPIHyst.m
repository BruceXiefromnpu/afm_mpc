

classdef TestPIHyst < matlab.unittest.TestCase
  
  properties
    
  end

  methods(Test)
    
    function test_invert_sat(self)
      t = (0:250)'*0.01;
      u = sin(t);
      
      ws = [0.0001; 0.02; .3; .5; 1.5; 0.03; .5; .065; 0.04];
      ymax = 1;
      
      n_d = 4;
      id_plus = (1:n_d);
      id_neg = (-n_d:-1);
      dplus = ((id_plus - 0.5)/n_d ) * ymax;
      dmin = ( (id_neg + 0.5)/n_d ) *ymax;
      d = [dmin, 0, dplus]';
      
      [dp, wsp] = PIHyst.invert_sat(d, ws);
      
      [~, usat] = PIHyst.sat_op(u, d, ws);
      
      [~, u_invsat_sat] = PIHyst.sat_op(usat, dp, wsp);
      
      self.verifyEqual(u, u_invsat_sat, 'AbsTol', 1e-9);
    end
    
    function test_invert_PI(self)
      nw = 7;
      umax = 2;
      r = ([0:nw-1]'./(nw) )*umax
      w = [0.0001; 0.02; .4; 0.78; 1; 0.03; 0.04];
      t = (0:250)'*0.01;
      u = sin(t);
      
      u_PI = PIHyst.hyst_play_op(u, r, w, w*0);
      [rp, wp] = PIHyst.invert_hyst_PI(r, w);
      u_invpi_pi = PIHyst.hyst_play_op(u_PI, rp, wp, wp*0);
      
      self.verifyEqual(u, u_invpi_pi, 'AbsTol', 1e-9)
    end
    
  end
  
  
end