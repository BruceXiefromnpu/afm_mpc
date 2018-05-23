classdef StageParams
%     This is a class of static parameters that define different aspects of our
%     stage and the sensors I have built, like the current measurement, and the
%     associated du_max derived from the TF bound.
% 
%     du_max
%     opamp2powI 
%     Vdiv2Vhigh 
%     Imax 
%     Rsense
  
  properties (Constant)
    du_max_1nm = 0.140;     % See publish_plotting/compute_du_1norm_bounds.m
    du_max_1nm_vib = 0.1007; % See publish_plotting/compute_du_1norm_bounds.m
    du_max = 0.198;
    opamp2powI = 1/15.15;
%     R1 = 29.7e6;
%     R2 = 1.732e6;
%     Vdiv_gain = R2/(R1 + R2);
    Vdiv2Vhigh = 1/(29.7e6/(1.732e6 + 29.7e6));
    Imax = 0.100;
    Rsense = 0.1;
    Ts = 40e-6;
  end

  methods
  end
  
end

