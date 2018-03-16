classdef condensedMPCprob_QP < CondensedMPCProb
%         H;
%         M;
%         Ainq;
%         binq;
%         N_mpc;
%         nu;              
%         ns;              
%         kappa;
%         warm_start_data;
%         sys;            
%
%
% Construction: condensedMPCprob(sys,N, Q,Qp, R)
%               condensedMPCprob(sys,N, Q,Qp, R, S)
%
% This sets H, M, N_mpc, and kappa fields. We leave it to you to fill in
% Ainq and binq via the method add_U_constraint().

    properties
        Ainq;            % Inequality LHS for input constraint.
        binq;            % Inequality RHS for input constraint. 
        lb = [];
        ub = [];
    end
    
    methods
        function self = condensedMPCprob_QP(sys,N, Q,Qp, R, S)
            % obj = condensedMPCprob(sys, N, Q, Qp, R, S)
            % Construct a conensedMPCprob instance. 
            if ~exist('S', 'var')
                ns = size(sys.b,1);
                nu = size(sys.b, 2);
                S = zeros(ns, nu);
            end
            [H, M] = CondensedMPCProb.build_mpc_problem(sys,N, Q, R, Qp, S);
            self.H     = H;
            self.M     = M;
            self.N_mpc = N;
            self.kappa = cond(H);
            self.nu = size(R,1);
            self.ns = size(Q,1);
            self.warm_start_data = [];
            self.sys = sys;
        end
        
        function U = call_qp_solver(self, f, Ainq, lbA, ubA, lb, ub)
            % lbA <= Ainq*U <= ubA
            %  is the same as
            %   [Ainq ] *U  <= [ubA ]
            %   [-Ainq]        [-lbA]
            Ainq = [self.Ainq;
                    Ainq;
                    -Ainq];
            binq = [self.binq;
                     ubA;
                     -lbA];
            opts.Display = 'off';  % Be quiet, quadprog.
            opts.MaxIter = 2000;
%             opts.ConstraintTolerance = 1e1;
            U = quadprog(self.H, f, Ainq, binq,[], [], lb, ub,...
                            [], opts);

        end
        
    end % METHODS
    
end % CLASSDEF






