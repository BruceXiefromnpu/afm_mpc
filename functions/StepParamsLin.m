classdef StepParamsLin 
    properties
        sys;
        ref_s;
        du_max;
        Q;
        S;
        gam_s;
        plant;
        trun;
        sim_struct;
    end
        
    methods
        function self = StepParamsLin(sys, ref_s, du_max, Q, gam_s, ...
                                      plant, trun, varargin)
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
            self.trun = trun;
            self.sim_struct = build_sim_struct(self);
            
        end % end contstructor
    end % end methods
end % end classdef
function sim_struct = build_sim_struct(self)
    xss = SSTools.getNxNu(self.plant);
    Nx = SSTools.getNxNu(self.sys);
    
    % Add a non-empty MPC problem, just so the simulink simulation
    % wont complain.
    mpcProb1 = condensedMPCprob(self.sys, 2, self.Q, self.Q, 100); 
    
    sim_struct = struct('PLANT', self.plant, 'trun', self.trun, 'mpcProb1', mpcProb1,...
                        'du_max', self.du_max, 'mpc_on', false,...
                        'xss', xss, 'Nx', Nx, 'N_mpc', 0);
    
    
end