classdef fict_zeros_test < matlab.unittest.TestCase
    
    properties
        sys;
    end
    
    methods 
       function self = fict_zeros_test()
            % Data from the power amp system.
            A = [1.429,   0.7208;
                -0.4429,   0.2865];
            B = [-0.2916;
                0.779];
            C = [0.3434,  0.1286];
            Ts = 40e-6;
            self.sys = ss(A, B, C, 0, Ts);
        end 
    end
    methods(Test)
        function test_fict_zeros(self)
        zdes = [0.5 + 0.5j, 0.5-0.5j];
        
        [Chat, Dhat, Q, S, Rmin] = fict_zeros(self.sys, zdes);
        
        sys_new = ss(self.sys.a, self.sys.b, Chat, Dhat, self.sys.Ts);
        
        z_act = tzero(sys_new)';
        z_act = sort(z_act);
        zdes = sort(zdes);

        self.verifyEqual(z_act, zdes, 'AbsTol', 1e-14);

        
        K = dlqr(self.sys.a, self.sys.b, Q, Rmin, S);
        
        a_cl = self.sys.a - self.sys.b*K;
        lqr_poles = sort(eig(a_cl)');
        
        self.verifyEqual(lqr_poles, zdes, 'AbsTol', 1e-14)
        
        
        end
    end
    
    
end