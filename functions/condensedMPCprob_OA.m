classdef condensedMPCprob_OA < handle
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
        H;                % The problem Hessian.
        M;                % The affine term, multiplied by xk_1.
%         Ainq;             % Inequality LHS for input constraint.
%         binq;             % Ineqauality RHS for input constraint.
%         lbAinq;           % Ineqaulity lower bound.
%         ubAinq;           % Inequality upper bound.
        CON;
        lb;               % Box constraint upper bound.
        ub;               % Box constraint lower bound.
        N_mpc;            % Control horizon.
        nu;               % Number of inputs.
        ns;               % Number of states.
        kappa;            % Condition number of H.
        warm_start_data;  % Solution from the previous iteration. This gets
                          % because we subclass the handle class.
        sys;              % The associated LTI ss system.
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
            [H, M] = clqrProblem_builder(sys,N, Q, R, Qp, S);
            obj.H     = H;
            obj.M     = M;
            obj.N_mpc = N;
            obj.kappa = cond(H);
            obj.nu = size(R,1);
            obj.ns = size(Q,1);
            obj.warm_start_data = qpOASES_auxInput('x0', zeros(N,1));
            obj.sys = sys;

        end
        
        function [U, X] = solve(self, xk_1, varargin)
        % Solve the condensed MPC problem using the qpOASES solver.
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
            
            if isempty(self.CON) 
                % || isempty(self.lbAinq) || isempty(self.ubAinq)
                [U,~,~,iter] = qpOASES(self.H, Mx0, self.lb, self.ub,{},...
                    self.warm_start_data);
            else  % Use the constraint class object. 
                self.CON.update_binq();
                
                [U,~,~,iter] = qpOASES(self.H, Mx0, self.CON.Ainq, self.CON.lb,...
                    self.CON.ub, self.CON.lbAinq, self.CON.ubAinq);
%                 , {},...
%                     self.warm_start_data);
                % Propogate the state forward.
                self.CON.update_sys(U(1:self.nu));
            end

            self.warm_start_data.x0 = [U(2:end);U(end)];
            
            if nargout == 2 
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
        % type: 'accum', accumlation constraint (useful for enforcing
        % magnitude saturation constraint when working with an
        % incremental form. bnds should be a vector.
            if strcmp('box', type) & length(bnds)==1
                bnds = [-bnds, bnds];
            end

            if strcmp(type, 'box')
                self.ub = ones(self.N_mpc, 1)*bnds(2);
                self.lb = ones(self.N_mpc, 1)*bnds(1);
            elseif strcmp(type, 'slew')
                S = derMat(self.N_mpc);
                self.Ainq = S;
                self.lbAinq = repmat(-bnds(1), self.N_mpc-1, 1);
                self.ubAinq = repmat(bnds(1), self.N_mpc-1, 1);
            elseif strcmp(type, 'accum')
               S = accumMat(self.N_mpc);
               self.Ainq = [self.Ainq; S];
               self.ubAinq = [self.ubAinq
                             repmat(bnds(2), self.N_mpc, 1)];
               self.lbAinq = [self.ubAinq
                             repmat(bnds(1), self.N_mpc, 1)];
                
            end
        end
        
%         function self = add_state_constraint(self, type, sys_aux, bnds)
%             
%             no = size(sys_aux.c, 1);
%             if strcmp('type', 'slew')
%                 f_0 = CondensedTools.init_cond_resp_matrix(sys_aux.a, 0,...
%                       self.N_mpc-no, sys_aux.c);
%                 f_1 = CondensedTools.init_cond_resp_matrix(sys_aux.a, 1,...
%                       self.N_mpc, sys_aux.c);
% 
%                 H_0 = CondensedTools.zero_state_output_resp(sys_aux,...
%                       self.N_mpc-1, 0);
% 
%                 H_1 = CondensedTools.zero_state_output_resp(sys_aux,...
%                       self.N_mpc, 1);
%                   
%                 self.Ainq = [self.Ainq; H_1-H_0];
%                 
%             else
%                 error('type not implmented')
%             end
%             
%         end
            
    end % METHODS

end % CLASSDEF
function M = accumMat(N)
    M = zeros(N, N);
    for k=1:N
       M(k, 1:k) = ones(1, k); 
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


function [H, M] = clqrProblem_builder(sys, N, Q, r, Qp, S)
    if sys.Ts == 0
        error('system "sys" should be discrete time dynamical system')
    end
    if isempty(Q) || isempty(r) || isempty(Qp) || isempty(S) ||isempty(N)
        error('System and weight matrices must not be empty!\n')
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
    %        A; 
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
    
    II_nplus1 = sparse(eye(N+1));
    II = sparse(eye(N));
    
    QQ = kron(II, Q);
    QQ = blkdiag(QQ, Qp);
    SS = kron(sparse(eye(N+1, N)), S);  % not the same as II_Nplus1 !!
    % RR = eye(N)*r;
    RR = kron(II, r);
    
    H = 2*(RR + F'*QQ*F + 2*F'*SS );

    H = (H+H')/2; %Symmetrize, should be anyway...
    M = 2*(phi'*QQ*F + phi'*SS)';


end

