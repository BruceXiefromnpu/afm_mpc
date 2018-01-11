classdef StepParamsTimeOpt 
    
    
    properties
        sys;
        ref_s;
        du_max;
        sys_nodelay;
        Nd;
    end
        
    methods
        function self = StepParamsTimeOpt(sys, ref_s, du_max, ...
                                          sys_nodelay, Nd)
            self.sys = sys;
            self.ref_s = ref_s;
            self.du_max = du_max;
            self.sys_nodelay = sys_nodelay;
            self.Nd = Nd;
        end % end contstructor
    end % end methods
end % end classdef
