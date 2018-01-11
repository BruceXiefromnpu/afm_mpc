classdef StepParamsCLQR 
    properties
        sys;
        ref_s;
        du_max;
        Q;
        gam_s;
        plant;
        N_traj;
        mpc_mode;
        
    end
        
    methods
        function self = StepParamsCLQR(sys, ref_s, du_max, Q, gam_s, ...
                                      plant, N_traj, mpc_mode)
            self.sys = sys;
            self.ref_s = ref_s;
            self.du_max = du_max;
            self.Q = Q;
            self.gam_s = gam_s;
            self.plant = plant;
            self.N_traj = N_traj;
            self.mpc_mode = mpc_mode;
        end % end contstructor
    end % end methods
end % end classdef
