classdef StepParamsMPC
    
    properties
        sys;
        ref_s;
        du_max;
        Q;
        R0;
        S;
        gam_s;
        plant;
        N_mpc_s;
        mpc_mode;
        trun;
        sim_struct;
        
    end
        
    methods
        function self = StepParamsMPC(sys, ref_s, du_max, Q, R0, gam_s, ...
                                      plant, trun, N_mpc_s, mpc_mode,...
                                      varargin)

            if ~strcmp(mpc_mode, 'sparse') && ~strcmp(mpc_mode, 'condensed')
                error('Expected mpc_mode = ("sparse"|"condensed"\n')
            end
            p = inputParser;
            p.addParameter('S', sys.b*0);
            p.parse(varargin{:});
            
            self.sys = sys;
            self.ref_s = ref_s;
            self.du_max = du_max;
            self.Q = Q;
            self.S = p.Results.S;
            self.R0 = R0;
            self.gam_s = gam_s;
            self.plant = plant;
            self.N_mpc_s = N_mpc_s;
            self.mpc_mode = mpc_mode;
            self.trun = trun;
            self.sim_struct = build_sim_struct(self);
            
        end
        

    end
    
end

function sim_struct = build_sim_struct(self)
    xss = SSTools.getNxNu(self.plant);
    Nx = SSTools.getNxNu(self.sys);
    
    sim_struct = struct('PLANT', self.plant, 'trun', self.trun, 'mpcProb1', [],...
                        'du_max', self.du_max, 'mpc_on', true,...
                        'xss', xss, 'Nx', Nx, 'N_mpc', 0);
    
    
end