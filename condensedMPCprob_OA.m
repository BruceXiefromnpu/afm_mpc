classdef condensedMPCprob_OA < handle
%         H;
%         M;
%         Ainq;
%         binq;
%         N_mpc;
%         kappa;
%
%
% Construction: condensedMPCprob(sys,N, Q,Qp, R)
%               condensedMPCprob(sys,N, Q,Qp, R, S)
%
% This sets H, M, N_mpc, and kappa fields. We leave it to you to fill in
% Ainq and binq.


    properties
        H;  % The problem Hessian.
        M;  % The affine term, multiplied by xk_1.
        Ainq; % Inequality LHS for input constraint.
        binq; % Ineqauality RHS for input constraint. 
        lb;
        ub;
        N_mpc; % Control horizon.
        nu;    % Number of inputs.
        ns;    % Number of states.
        kappa; % Condition number of H.
        n_warmstart; % set to N_mpc right now.
        warm_start_data;
        sys;         % The associated LTI ss system.
%         n_dec_vec;
    end
    
    methods
        function obj = condensedMPCprob_OA(sys,N, Q,Qp, R, S)
            % obj = condensedMPCprob(sys,N, Q,Qp, R, S)
            % Construct a conensedMPCprob instance. 
            if isa(sys, 'condensedMPCprob')
                obj.H = sys.H;
                obj.M = sys.M;
                obj.sys = sys;
                obj.kappa = sys.kappa;
                obj.Ainq = sys.Ainq;
                obj.binq = sys.binq;
                obj.N_mpc = sys.N_mpc;
                obj.nu = sys.nu;
                obj.ns = sys.ns;
                obj.n_warmstart = sys.n_warmstart;
                obj.warm_start_data = sys.warm_start_data;
                obj.n_dec_vec = sys.n_dec_vec;
            else
                if ~exist('S', 'var')
                    ns = size(sys.b,1);
                    nu = size(sys.b, 2);
                    S = zeros(ns, nu);
                end
                [H, M] = clqrProblem_local(sys,N, Q, R, Qp, S);
                obj.H     = H;
                obj.M     = M;
                obj.N_mpc = N;
                obj.kappa = cond(H);
                obj.nu = size(R,1);
                obj.ns = size(Q,1);
                obj.n_warmstart = N;
                obj.warm_start_data = qpOASES_auxInput('x0', zeros(N,1));
                obj.sys = sys;
                
            end
        end
        
        function [U, X] = solve(self, xk_1, varargin)
            % Solve the condensed MPC problem.
            % [U, X] = solve(self, xk_1, varargin)
            %
            % inputs
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
            p.addParameter('getX', false);
            parse(p, varargin{:});
            getX = p.Results.getX;
            Mx0 = self.M*xk_1;

            [U,~,~,iter] = qpOASES(self.H, Mx0, self.lb, self.ub,{},...
                self.warm_start_data);

            self.warm_start_data.x0 = [U(2:end);U(end)];
            
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
           
            % uk = Us(1:self.nu);
        end
        
        function self = add_U_constraint(self, type, bnds)
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
                self.ub = ones(self.N_mpc, 1)*bnds(2);
                self.lb = ones(self.N_mpc, 1)*bnds(1);
%                 I = eye(self.N_mpc*self.nu);
%                 Ainq = [I;-I];
%                 binq = [ones(self.N_mpc,1)*bnds(2);
%                        ones(self.N_mpc,1)*(-bnds(1))];
            elseif strcmp(type, 'slew')
                S = derMat(self.N_mpc);
                self.Ainq = [S; -S];
                self.binq = [zeros(2*self.N_mpc-2, 1)+bnds(1)];
            end

        end
        
        
    end
    
end

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
    if sys.Ts == 0
        error('system "sys" should be discrete time dynamical system')
    end

    Ns = size(sys.b, 1);
    nu = size(sys.b, 2);

    [a, b, c, d] = ssdata(sys);

    % -------------------------------------------------------------------------
    % Condense the LQR probelm
    % X(k) = phi*x(0) + F*[u(0)
    %                      u(1)
    %                       :
    %                       u(N-1)]
    % Where G(i,:) = C*A^i
    % and 


    % make phi;
    % phi = [
%            A; 
    %       A^2;
    %       A^3;
    %         :
    %       A^N]

    phi = []; 
    F = [];
    for i = 1:N
        phi = [phi; a^i];
    end

    F_1 = [eye(Ns); phi(1:end-Ns, :)]*b;

    for i=0:N-1
       F = [F, [zeros(Ns*i,nu); F_1(1:Ns*N-(i)*Ns, :)] ];
    end

    phi = [eye(Ns);
           phi];
    F = [zeros(Ns, nu*N); F];
    
    
    
    II_nplus1 = eye(N+1);
    QQ = kron(II_nplus1, Q);
    QQ = kron(eye(N), Q);
    QQ = blkdiag(QQ, Qp);
%     QQ(N*Ns+1-Ns:end, N*Ns+1-Ns:end) = Qp;

    SS = kron(eye(N+1, N), S);

    % RR = eye(N)*r;
    RR = kron(eye(N), r);
    
    H = 2*(RR + F'*QQ*F + 2*F'*SS );

    H = (H+H')/2; %Symmetrize, should be anyway...
    M = 2*(phi'*QQ*F + phi'*SS)';


end

