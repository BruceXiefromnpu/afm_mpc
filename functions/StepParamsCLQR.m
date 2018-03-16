classdef StepParamsCLQR 
    properties
        sys;
        ref_s;
        du_max;
        Q;
        S;
        gam_s;
        plant;
        N_traj;
        mpc_mode;
        
    end
        
    methods
        function self = StepParamsCLQR(sys, ref_s, du_max, Q, gam_s, ...
                                      plant, N_traj, mpc_mode, varargin)

            if ~strcmp(mpc_mode, 'sparse') && ~strcmp(mpc_mode, 'condensed')
                error('Expected mpc_mode = ("sparse"|"condensed"\n')
            end
            p = inputParser;
            p.addParameter('S', sys.b*0);
            p.parse;
            self.sys = sys;
            self.ref_s = ref_s;
            self.du_max = du_max;
            self.Q = Q;
            self.S = p.Results.S;
            self.gam_s = gam_s;
            self.plant = plant;
            self.N_traj = N_traj;
            self.mpc_mode = mpc_mode;
        end % end contstructor
    end % end methods
end % end classdef
