classdef sparseMPCprob < handle
%   sparseMPCprob  container for the sparse MPC QP formulation.
%         H;      The problem Hessian
%         Aeq;    The LHS equality constraint enforcing dynamcs.
%         beq;    the RHS equality equality constraint enforcing dynamcs.
%         Ainq;   Inequality LHS for control constraint.  
%         binq;   Inequality RHS for control constraint.  
%         N_mpc;  Control horizon
%         kappa;  condition number of H
%         ns;     number of states.
%         nu;     Number of controls.
%         qpopts; Options that get passed to quadrog. From optimset.
    %
    % Construction: condensedMPCprob(sys,N, Q,Qp, R)
    %               condensedMPCprob(sys,N, Q,Qp, R, S)
    %
    % This sets H, M, N_mpc, and kappa fields. We leave it to you to fill in
    % Ainq and binq.


    properties
        H;      % The problem Hessian
        Aeq;    % The LHS equality constraint enforcing dynamcs.
        beq;    % the RHS equality equality constraint enforcing dynamcs.
        Ainq;   % Inequality LHS for control constraint.  
        binq;   % Inequality RHS for control constraint.  
        N_mpc;  % Control horizon
        kappa;  % condition number of H
        ns;     % number of states.
        nu;     % Number of controls.
        qpopts; % Options that get passed to quadrog. From optimset.
    end
    
    methods
        function obj = sparseMPCprob(sys,N, Q,Qp, R, S)
            % obj = sparseMPCprob(sys,N, Q,Qp, R, S)
            % S can be [], in which case, it is set to the zero matrix.
            if ~exist('S', 'var')
                ns = size(sys.b,1);
                nu = size(sys.b, 2);
                S = zeros(ns, nu);
            end
            [H, Aeq, beq] = clqrProblem_local(sys,N, Q, R, Qp, S);
            obj.H     = sparse(H);
            obj.Aeq   = sparse(Aeq);
            obj.beq   = beq;
            obj.Ainq  = [];
            obj.binq  = [];
            obj.N_mpc = N;
            obj.kappa = cond(H);
            obj.ns = size(sys.b,1);
            obj.nu = size(sys.b,2);
            obj.qpopts = optimset('Display', 'off');
        end
        
        
        function [U, X] = solve(self, xk_1, varargin)
            % [U, X] = solve(self, xk_1)
            % Solves the MPC problem defined by the instance, given
            % an initial condition xk_1. 
            % 
            % Inputs
            % ------
            %   xk_1 : the initial condition at which to solve the MPC
            %   problem
            %
            %   solve(..., 'getX', {0|1}) If 1, also return the state
            %   sequence. Since we are in the condensed version, this is
            %   found via lsim(...)
            %  
            %   solve(...,'warm_start_data', not used currently.
            %
            % Outputs
            % ------
            %   U : the sequence of optimal controls from k=1...self.N_mpc
            %   X : the sequence of optimal states from k= 1...self.N_mpc.
            %        Empty if 'getX' is 0 (default), otherwise, a matrix 
            %        of dimensions (Ns x N_mpc).
            
            % Parse varagin.
            defaultGetX = true;
            defaultWarm_start_data = [];
            p = inputParser;
            p.addParameter('getX', defaultGetX);
            p.addParameter('warm_start_data', defaultWarm_start_data);
            parse(p, varargin{:});
            getX = p.Results.getX;
            
            beq_xk = self.beq;
            beq_xk(1:self.ns) = xk_1;
            UX = mpcSolve_local(self.H, self.Ainq, self.binq,...
                self.Aeq, beq_xk, self.qpopts);
            
            % Split the solution into controls and states.
            U = UX(1:self.nu*self.N_mpc)';
            if getX
                X = UX(self.nu*self.N_mpc+1:end);
                X = reshape(X, self.ns, []);
            else
                X = [];
            end
        end
        
        function self = add_U_constraint(self, type, bnds)
            % self = add_U_constraint(self, type, bnds)
            %
            % Add an input constraint to the MpcPRoblem instance.
            % type: 'box' forms a box constraint on u between bnds(1) and
            % bnds(2)
            % type: symmetric 'slew' slew rate constraint. bnds is scalar
            if strcmp('box', type) & length(bnds)==1
                bnds = [-bnds, bnds];
                fprintf('Type box constraint indicated but only one bound found\n')
                fprintf('Assuming symmentric bounds of [%.3f, %.3f]\n', bnds(1), bnds(2));
            end

            if strcmp(type, 'box')
                Zro = zeros(self.N_mpc*self.nu*2, (self.N_mpc+1)*self.ns);
                I = eye(self.N_mpc*self.nu);
                Ainq = [[I;-I], Zro];
                binq = [ones(self.N_mpc,1)*bnds(2);
                       ones(self.N_mpc,1)*(-bnds(1))];
            elseif strcmp(type, 'slew')
                Zro = zeros((self.N_mpc*self.nu-1)*2, (self.N_mpc+1)*self.ns);
                S = derMat(self.N_mpc);
                Ainq = [[S; -S], Zro];
                binq = [zeros(2*self.N_mpc-2, 1)+bnds(1)];
            end
            self.Ainq = [self.Ainq; Ainq];
            self.binq = [self.binq; binq];
        end
        
    end % METHODS
end % CLASDEF

function UX = mpcSolve_local(H, Ainq,binq, Aeq, beq, qpopts)

    f  = zeros(size(H,2),1);
    UX = quadprog(H, f, Ainq, binq, Aeq, beq, [], [], [], qpopts);
end


function [H, Aeq, beq] = clqrProblem_local(sys, N, Q, R, Qp, S)
    if sys.Ts == 0
        error('system "sys" should be discrete time dynamical system')
    end
    
    ns = size(sys.b, 1);
    nu = size(sys.b, 2);

    [a, b] = ssdata(sys);

    % -------------------------------------------------------------------------
    % Sparse cLQR problem

    % Build the cost function. THis should look like:
    %      [R    |S 0 0 ][u0]
    %      [   R |0 S 0 ][u1]
    %  H = [-----|------][--]
    %      [S' 0 |Q     ][x0]
    %      [0  S'|  Q   ][x1]
    %      [0  0 |    Qp][xN]

    % First, build the part without the cross term. Use eye(N) because
    % 0...N-1 has N elements.

    RRR = kron(eye(N), R);
    QQQp = blkdiag(kron(eye(N), Q), Qp);
    Sx = [kron(eye(N), S);
      zeros(ns, N*nu)]; % The lower left.
    Su = Sx';             % The upper right.

    % The main event:
    H = [RRR, Su;
        Sx, QQQp];
    
    
    % Now build the equality constraint. This should look like:
    %  [I 0 0  0 0 0 ][x0]   [x0]
    %  [A B -I 0 0 0 ][u0] = [0 ]
    %  [0 0  A B -I  ][x1]   [0 ]
    %                 [u1]
    %                 [x_N] 
    
    
    I_b = eye(N, N);
    Aeq_b = kron(I_b, b);
    Aeq_a = [kron(eye(N), a), zeros(N*ns, ns)]...
            + [zeros(N*ns, ns), kron(eye(N), -eye(ns))];
    
    
    Aeq = [[zeros(ns, nu*N), eye(ns), zeros(ns, N*ns)];
           [Aeq_b, Aeq_a]];
           

    beq = zeros(ns*(N+1), 1);

end


function S = derMat(N)
% S = derMat(N)
% Computes and N-1 x N derivitave matrix.

    S = -eye(N);
    I2 = eye(N-1);

    S(1:N-1, 2:N) = I2+S(1:N-1, 2:N);

    S = S(1:N-1, :);
end

