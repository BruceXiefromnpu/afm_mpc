classdef CanonRefTraj
  
  methods(Static)

    function [ref_traj, trun] = ref_traj_1(r1, N)
      N1 = N;
      N2 = N;
      Ts = StageParams.Ts;
      t1 = [0:1:N1]'*Ts;
      trun = Ts*(N1);
      t = t1;
      yref = [0*t1 + r1];
      ref_0 = 0;
      ref_traj = timeseries(yref, t);      
      
    end
    
    function [ref_traj, trun] = ref_traj_2(r1, r2, N)
      N2 = N;
      N1 = N;
      
      trun = StageParams.Ts*(N1 + N2);
      ref_f_2 = 1.5; % 1.5 to hit slew rate, 1.4 doesn't
      ref_0 = 0;
      t1 = [0:1:N1]'*StageParams.Ts;
      t2 = [N1+1:1:N1+N2]'*StageParams.Ts;
      
      t = [t1; t2];
      yref = [0*t1 + r1;
              0*t2 + r2];
      ref_traj = timeseries(yref, t);      
      
    end
    
    function ref_traj = ref_traj_load(fpath)
      load(fpath)
      ref_traj = ref_traj_params.ref_traj;
    end
    
  end
  
end

  
  
