classdef condensedMPCprob < handle
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
        H;               % The problem Hessian.
        M;               % The affine term, multiplied by xk_1.
        Ainq;            % Inequality LHS for input constraint.
        binq;            % Inequality RHS for input constraint. 
        N_mpc;           % Control horizon.
        nu;              % Number of inputs.
        ns;              % Number of states.
        kappa;           % Condition number of H.
        warm_start_data; % Data from prior iteration for warmstarting next.
        sys;             % The associated LTI ss system.
    end
    
    methods
        function self = condensedMPCprob(sys,N, Q,Qp, R, S)
            % obj = condensedMPCprob(sys, N, Q, Qp, R, S)
            % Construct a conensedMPCprob instance. 
            if ~exist('S', 'var')
                ns = size(sys.b,1);
                nu = size(sys.b, 2);
                S = zeros(ns, nu);
            end
            [H, M] = clqrProblem_local(sys,N, Q, R, Qp, S);
            self.H     = H;
            self.M     = M;
            self.N_mpc = N;
            self.kappa = cond(H);
            self.nu = size(R,1);
            self.ns = size(Q,1);
            self.warm_start_data = [];
            self.sys = sys;
        end
        
        function [U, X] = solve(self, xk_1, varargin)
            % Solve the condensed MPC problem.
            % [U, X] = solve(self, xk_1, varargin)
            %
            % Inputs
            % ------
            %   xk_1 : current state (ie, initial condition) for the
            %   optimization problem.
            %
            %   solve(..., 'getX', {0|1}) If 1, also return the state
            %   sequence. Since we are in the condensed version, this is
            %   found via lsim(...)
            %  
            %   solve(...,'warm_start_data', not used currently.
            %
            % Outputs
            % ------
            %   U : vector of optimal controls
            %   X : Empty if getX is 0 (default), otherwise, a matrix of
            %   optmial states of size ns x N_mpc.
            
            % Parse varagin.
            p = inputParser;
            p.addParameter('getX', 0);
            
            parse(p, varargin{:});
            getX = p.Results.getX;
            
            f = self.M*xk_1;
            
            opts.Display = 'off';  % Be quiet, quadprog.
            U = quadprog(self.H, f, self.Ainq, self.binq,[], [], [], [],...
                            self.warm_start_data, opts);
            self.warm_start_data = [U(2:end); U(end)];
            
            if getX 
               X = discrete_lsim(self.sys, U,  xk_1);
            else
               X = [];
            end
            
            % For consistency with the shape of X, and with sparseMPCprob,
            % reshape so that cols correspond to time k, and rows are input
            % channel
            % [ u_1(0) u_1(0) ...]
            % [ u_2(0) u_2(1) ...]
            U = reshape(U, self.nu, []);
        end
        
        function self = add_U_constraint(self, type, bnds)
            % Add an input constraint to the MpcPRoblem instance.
            % type: 'box' forms a box constraint on u between bnds(1) and
            % bnds(2)
            % type: symmetric 'slew' slew rate constraint. bnds is scalar
            if strcmp('box', type) && length(bnds)==1
                bnds = [-bnds, bnds];
                fprintf('Type box constraint indicated but only one bound found\n')
                fprintf('Assuming symmentric bounds of [%.3f, %.3f]\n', bnds(1), bnds(2));
            end
            
            if strcmp(type, 'box')
                I = eye(self.N_mpc*self.nu);
                Ainq = [I;-I];                             %#ok<PROPLC>
                binq = [ones(self.N_mpc,1)*bnds(2);        %#ok<PROPLC>
                       ones(self.N_mpc,1)*(-bnds(1))];
            elseif strcmp(type, 'slew')
                S = derMat(self.N_mpc);
                Ainq = [S; -S];                            %#ok<PROPLC>
                binq = [zeros(2*self.N_mpc-2, 1)+bnds(1)]; %#ok<PROPLC>
            end
            self.Ainq = [self.Ainq; Ainq];                 %#ok<PROPLC>
            self.binq = [self.binq; binq];                 %#ok<PROPLC>
        end
    end % METHODS
    
end % CLASSDEF


function X = discrete_lsim(sys, U, x0)
    % The sparseMPC returns the N+1th state. This function does the same 
    % to maintain consistency between the two classes.
    xk = x0;
    X = zeros(size(sys.B, 1), size(U,1)+1);
    for k = 1:length(U)
        uk = U(k, :);
        X(:,k) = xk;
        xk_1 = sys.A*xk + sys.B*uk;
        xk = xk_1;
    end
    X(:,k+1) = xk_1;

end


function [H, M] = clqrProblem_local(sys, N, Q, r, Qp, S)
% Condensed MPC problem
% X(k) = phi*x(0) + F*[u(0)
%                      u(1)
%                       :
%                       u(N-1)]
% Where F(i,:) = A^iB
% and 
% phi = [
%        A; 
%       A^2;
%       A^3;
%         :
%       A^N]
    
    if sys.Ts == 0
        error('system "sys" should be discrete time dynamical system')
    end
    
    ns = size(sys.b, 1);
    nu = size(sys.b, 2);

    [a, b] = ssdata(sys);

    phi = [];  % Should pre-allocate but I was lazy. 
    F = [];
    for i = 1:N
        phi = [phi; a^i];
    end

    F_1 = [eye(ns); phi(1:end-ns, :)]*b;

    for i=0:N-1
       F = [F, [zeros(ns*i,nu); F_1(1:ns*N-(i)*ns, :)] ];
    end

    phi = [eye(ns);
           phi];
    F = [zeros(ns, nu*N); F];
    
    II = sparse(eye(N));
    II_Nplus1 = sparse(eye(N+1, N));
    
    QQ = kron(eye(N), Q);
    QQ = blkdiag(QQ, Qp);
    SS = kron(II_Nplus1, S);
    RR = kron(II, r); % if r scalar, RR = I*r.
    
    H = 2*(RR + F'*QQ*F + 2*F'*SS );

    H = (H+H')/2; %Symmetrize, should be anyway...
    M = 2*(phi'*QQ*F + phi'*SS)';


end

