classdef CondensedMPCProb < handle
    %CONDENSEDMPC Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        H;                % The problem Hessian.
        M;                % The affine term, f = M*xk_1.
        CON;              % A constraint object.
        N_mpc;            % Control horizon.
        nu;               % Number of inputs.
        ns;               % Number of states.
        kappa;            % Condition number of H.
        warm_start_data;  % Solution from the previous iteration. This gets
                          % because we subclass the handle class.
                          
        sys;              % The associated LTI ss system.
    end
    
    methods (Abstract)
        % Solve the mpc problem for initial condition xk_1.
        call_qp_solver(self, f, Ainq, lbA, ubA)
        
        %  add_U_constraint(self, type, bnds)
        % Add an input constraint to the MpcPRoblem instance.
        % type: 'box' forms a box constraint on u between bnds(1) and
        % bnds(2)
        % type: symmetric 'slew' slew rate constraint. bnds is scalar
        % type: 'accum', accumlation constraint (useful for enforcing
        % magnitude saturation constraint when working with an
        % incremental form. bnds should be a vector.
    end
    
    methods
        function [U, X] = solve(self, xk_1, varargin)
        % [U, X] = solve(self, xk_1, varargin)
        % Solve the mpc problem for initial condition xk_1.
            f = self.M*xk_1;
            
            if isempty(self.CON) 
                Ainq = [];
                lbA = [];
                ubA = [];
                [U] = self.call_qp_solver(f, Ainq, lbA, ubA, self.lb, self.ub);
            else  % Use the constraint class object. 
                lb = self.CON.lb;
                ub = self.CON.ub;

                self.CON.update_binq();

                U = self.call_qp_solver(f, self.CON.Ainq,...
                                        self.CON.lbAinq,...
                                        self.CON.ubAinq, lb, ub);
                
                self.CON.update_sys(U(1:self.nu));
            end
            % For consistency with the shape of X, and with sparseMPCprob,
            % reshape so that cols correspond to time k, and rows are input
            % channel
            % [ u_1(0) u_1(0) ...]
            % [ u_2(0) u_2(1) ...]
            U = reshape(U, self.nu, []);
        
            if nargout == 2
               X = self.discrete_lsim(U',  xk_1);
            else
               X = [];
            end
        end
% %         function self = add_U_constraint(self, type, bnds)
% %         % Add an input constraint to the MpcPRoblem instance.
% %         % type: 'box' forms a box constraint on u between bnds(1) and
% %         % bnds(2)
% %         % type: symmetric 'slew' slew rate constraint. bnds is scalar
% %         % type: 'accum', accumlation constraint (useful for enforcing
% %         % magnitude saturation constraint when working with an
% %         % incremental form. bnds should be a vector.
% %             if strcmp('box', type) & length(bnds)==1
% %                 bnds = [-bnds, bnds];
% %             end
% % 
% %             if strcmp(type, 'box')
% %                 self.ub = ones(self.N_mpc, 1)*bnds(2);
% %                 self.lb = ones(self.N_mpc, 1)*bnds(1);
% %             elseif strcmp(type, 'slew')
% %                 S = derMat(self.N_mpc);
% %                 self.Ainq = S;
% %                 self.lbAinq = repmat(-bnds(1), self.N_mpc-1, 1);
% %                 self.ubAinq = repmat(bnds(1), self.N_mpc-1, 1);
% %             elseif strcmp(type, 'accum')
% %                S = self.accumMat(self.N_mpc);
% %                self.Ainq = [self.Ainq; S];
% %                self.ubAinq = [self.ubAinq
% %                              repmat(bnds(2), self.N_mpc, 1)];
% %                self.lbAinq = [self.ubAinq
% %                              repmat(bnds(1), self.N_mpc, 1)];
% %                 
% %             end
% %         end        
        function X = discrete_lsim(self, U, x0)
        % The sparseMPC returns the N+1th state. This function does the same 
        % to maintain consistency between the two classes.
            xk = x0;
            X = zeros(size(self.sys.B, 1), size(U,1)+1);
            for k = 1:length(U)
                uk = U(k, :);
                X(:,k) = xk;
                xk_1 = self.sys.A*xk +self.sys.B*uk;
                xk = xk_1;
            end
            X(:,k+1) = xk_1;
        end % end discrete_lsim        

        function M = accumMat(~, N)
        % Builds an accumulator matrix M such that for a vector U, M*U is
        % equivalent to cumsum(U).
            M = zeros(N, N);
            for k=1:N
                M(k, 1:k) = ones(1, k);
            end
        end
        
    end % End methods
    
    methods (Static)
        function [H, M] = build_mpc_problem(sys, N, Q, r, Qp, S)
        % Condense the LQR problem
        % | x0 | = | I  |    |0  0  0 0 0 | | u(0) |
        % | x1 |   | A  |    |B  0  0 0 0 | | u(2) |
        % |  : |   |A^2 |x0 +|AB B  0 0 0 | | u(3) |
        % |  : |   |A^3 |    |A^2B AB B 0 | |  :   |
        % |    |   | :  |    |A^N-1B ...B | |u(N-1)|
        % | xN |   |A^N |      
        %   ^        ^             ^          ^
        %   |        |             |          |
        %   X       AA             BB         U
            
            if sys.Ts == 0
                error('system "sys" should be discrete time dynamical system')
            end
            if isempty(Q) || isempty(r) || isempty(Qp) || isempty(S) ||isempty(N)
                error('System and weight matrices must not be empty!\n')
            end
            Ns = size(sys.b, 1);
            nu = size(sys.b, 2);

            [a, b, c, d] = ssdata(sys);

            AA = []; 
            BB = [];
            for i = 0:N
                AA = [AA; a^i];
            end
            
            % The first column of BB. Dont want the A^N term.
            BB_col1 = AA(1:end-Ns, :)*b;            
            % Build column-by-column
            for i=0:N-1
               BB = [BB, [zeros(Ns*i,nu); BB_col1(1:Ns*N-(i)*Ns, :)] ];
            end
            % Pad with zeros at the top.
            BB = [zeros(Ns, nu*N); BB];

            II = sparse(eye(N));

            QQ = kron(II, Q);
            QQ = blkdiag(QQ, Qp);
            SS = kron(sparse(eye(N+1, N)), S);  % not the same as II_Nplus1 !!
            RR = kron(II, r);  % RR = eye(N)*r;
            
            % N.B.  BB'*SS + SS'*BB !=2BB'*SS
            H = 2*(RR + BB'*QQ*BB + BB'*SS + SS'*BB); 
            H = (H+H')/2; %Symmetrize, should be anyway...

            M = 2*(AA'*QQ*BB + AA'*SS)';


        end % end builder
    end % end methods(Static)
end

