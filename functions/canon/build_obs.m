function [sys_obsDist, L_dist] = build_obs(sys_obs, can_params)
 %  [sys_obsDist, L_dist] = build_obs(sys_obs, can_params)
  if sys_obs.InputDelay ~= 0 || sys_obs.IOdelay ~=0
    error(['All Delays of sys_obs must be part of the states, not ' ...
           'included as a property']);
  end
  
  
  p_int_d = can_params.p_int;
  
  Qw = sys_obs.b*sys_obs.b'*can_params.beta;
  
  Lx = dlqr(sys_obs.a', sys_obs.c', Qw, 1)';
  [L_dist, sys_obsDist, IDENT_obs, eNs_12] = DistEst.output_dist_est(sys_obs,...
                                                    Lx, p_int_d);
  [Nx_r, Nx_d, Nu_r, Nu_d] = DistEst.steady_state_gains(sys_obs, sys_obs.b*0, 1);
  
end
