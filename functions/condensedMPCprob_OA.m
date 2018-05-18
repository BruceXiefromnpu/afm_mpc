classdef condensedMPCprob_OA < CondensedMPCProb
%         H;                The problem Hessian.
%         M;                The affine term, multiplied by xk_1.
%         Ainq;             Inequality LHS for input constraint.
%         binq;             Ineqauality RHS for input constraint. 
%         lb;               Box constraint upper bound.
%         ub;               Box constraint lower bound.
%         N_mpc;            Control horizon.
%         nu;               Number of inputs.
%         ns;               Number of states.
%         kappa;            Condition number of H.
%         warm_start_data;  Solution from the previous iteration. This gets
%                           because we subclass the handle class.
%         sys;              The associated LTI ss system.
%
%
% Construction: condensedMPCprob(sys,N, Q,Qp, R)
%               condensedMPCprob(sys,N, Q,Qp, R, S)
%
% This sets H, M, N_mpc, and kappa fields. We leave it to you to fill in
% Ainq and binq.


    properties
%         H;                % The problem Hessian.
%         M;                % The affine term, multiplied by xk_1.
%         CON;
        lb;               % Box constraint upper bound.
        ub;               % Box constraint lower bound.
%         N_mpc;            % Control horizon.
%         nu;               % Number of inputs.
%         ns;               % Number of states.
%         kappa;            % Condition number of H.
%         warm_start_data;  % Solution from the previous iteration. This gets
%                           % because we subclass the handle class.
%         sys;              % The associated LTI ss system.
    end
    
    methods
        function obj = condensedMPCprob_OA(sys,N, Q,Qp, R, S)
            % obj = condensedMPCprob(sys,N, Q,Qp, R, S)
            % Construct a condensedMPCprob_OA instance. 
            if ~exist('S', 'var')
                ns = size(sys.b,1);
                nu = size(sys.b, 2);
                S = zeros(ns, nu);
            end
            [H, M] = CondensedMPCProb.build_mpc_problem(sys,N, Q, R, Qp, S);
            obj.H     = H;
            obj.M     = M;
            obj.N_mpc = N;
            obj.kappa = cond(H);
            obj.nu = size(R,1);
            obj.ns = size(Q,1);
            obj.warm_start_data = qpOASES_auxInput('x0', zeros(N,1));
            obj.sys = sys;
            obj.ub = []; 
            obj.lb = []; 

        end
        function reset_warm_start_data(self)
          N = self.N_mpc;
          self.warm_start_data = qpOASES_auxInput('x0', zeros(N,1));
        end
        function U = call_qp_solver(self, f, Ainq, lbA, ubA, lb, ub)
            
            if isempty(Ainq) 
                % || isempty(self.lbAinq) || isempty(self.ubAinq)
                [U,~,exit_flag] = qpOASES(self.H, f, lb, ub,{},...
                    self.warm_start_data);
            else  % Use the constraint class object. 

                [U,~,exit_flag] = qpOASES(self.H, f, Ainq, lb,...
                    ub, lbA, ubA); %, {}, self.warm_start_data);
                if exit_flag ~=0
                    keyboard
                   condensedMPCprob_OA.handle_exitFlag(exit_flag) 
                end
            end

            self.warm_start_data.x0 = [U(2:end);U(end)];
            
        end
        
    end % METHODS
    
    methods (Static)
   
        function handle_exitFlag(exit_flag)
           if exit_flag == 1
               fprintf('QP could not be solved within given number of iterations\n');
           elseif exit_flag == -1
               fprintf('QP could not be solved due to an internal error\n');
           elseif exit_flag == -2
               fprintf('QP is infeasible (and thus could not be solved)\n');
           elseif exit_flag == -3
               fprintf('QP is unbounded (and thus could not be solved)\n');
           end
        end
    end % methods(Static)

end % CLASSDEF






