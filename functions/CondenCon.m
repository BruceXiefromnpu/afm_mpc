classdef CondenCon < handle
% This class provides an abstraction layer for more generalized constraint
% handling for optimizations that employ the condensed format. The goal is
% that this class should be usable by both condensedMPCprob based classes as
% well as the TimeOptBisect class (which also optimizes over a vector of U,
% and needs constraints that are of the same for as for condensedMPCprob).
%
% The difficulty with state constraints is that vector sides of the linear
% inequality depend on the initial condition of the state at time k 
%            lb(x_k) <= Ainq * U <= ub(x_k)
%
% For now, I don't want to break the API that I have developed for my 
% single simulink model setup. That may have to happen in the future if I
% wanted to use a state *estimate* of the constrained state. For now, I
% will propogate that state in open loop inside this class. This shouldn't
% matter for my initial simulations and its unclear to me right now exactly
% how that state estimation problem would work if I don't either (1) add a
% sensor measurement at the power amp output or (2) force my overall system
% model to explicitly include the power amp model.
%
% Construction Inputs
% -------------------
%   sys : a state space system representing the state which should be
%         constrained.
%
%   x0 :  The initial condition of sys.
%
%   N_mpc : the control horizon.
%
% Usage
% -----
%   The primary use-case is for state constraints. I illustrate that first.
% 
%   1. Construction 
%   >> con = CondenCon(sys, x0, N)
%   where sys is the state space system which contains the states to
%   constrain. The constraint will be over the ouput y of sys.
%
%  2. Add constraints.
%
%
%
% TODO
% ----
% ** Support more general state constraints that a SISO output y. The main
%    change that is needed here is handled the bnds properly in
%    add_state_con().
% 
% ** Fully support the 'accum' option. Basically, we need to accumulate u
%    inside this class, and the update the constraint based on the previous
%    value of u.

    properties
        lbAinq_funs;  % A cell array of function handles that update the LHS
                      % of the inequality.
        ubAinq_funs;  % A cell array of function handles that update the RHS
                      % of the inequality.
        lbAinq;       % The current value of the LHS inequality.
        ubAinq;       % The current value of the RHS inequality.
        lb;           % Lower bound constraint on input.
        ub;           % Upper bound constraint on input.
        Ainq;         % The linear constraint matrix in 
        sys;
        x0;
        N_mpc
        xvec;
    end
    
    methods
        function self = CondenCon(sys, x0, N_mpc)
            self.sys = sys;
            self.x0 = x0;
            self.xvec = x0;
            self.N_mpc = N_mpc;
            self.lbAinq_funs = {};
            self.ubAinq_funs = {};
        end
        
        function update_sys(self, uk)
            if ~isempty(self.sys)
               self.x0 = self.sys.a*self.x0 + self.sys.b*uk;
%                xk = self.sys.a*self.xvec(:,end) + self.sys.b*uk;
               self.xvec = [self.xvec, self.x0];
            end
        end
        
        function update_binq(self)
        % update_binq()
        % 
        % Updates the upper and lower bound of the linear inequality by
        % evaluating the array of function handles in
        % self.(lb)(ub)Ainq_funs.
        %
            self.ubAinq = [];
            self.lbAinq = [];
            for k=1:length(self.lbAinq_funs)
                self.lbAinq = [self.lbAinq; 
                    self.lbAinq_funs{k}(self.x0)];
            end
            for k=1:length(self.ubAinq_funs)
                self.ubAinq = [self.ubAinq; 
                    self.ubAinq_funs{k}(self.x0)];
            end
            
        end
        
        function self = add_state_con(self, type, bnds)
        % self = add_state_con(self, type, bnds)
        % Adds a state constraint to a CondenCon instance. 
        % 
        % Inputs
        % ------
        % type : string, 'box' is the only currently supported type.
        % 
        % bnds : is a 1x2 vector of upper and lower bounds. 
            no = size(self.sys.c, 1);
            if strcmp(type, 'box')
                f_1 = CondenCon.init_cond_resp_matrix(self.sys.a, 0,...
                      self.N_mpc, self.sys.c);

                H_1 = CondenCon.zero_state_output_resp(self.sys,...
                      self.N_mpc, 0);
                H_1 = H_1(:,1:end-1);  
                self.Ainq = [self.Ainq; H_1];
                ymax_s = ones(size(f_1,1), 1)*bnds(1);
                self.ubAinq_funs{end+1} = @(x0) ymax_s - f_1*x0;
                self.lbAinq_funs{end+1} = @(x0) -ymax_s - f_1*x0;
            else
                error('type not implmented')
            end        
        end
        
        % add_state_con(self, type, bnds)
        function add_input_con(self, type, bnds)
        % add_state_con(self, type, bnds)
        %
        % Adds an input constraint to a CondenCon instance.
        %
        % Inputs
        % ------
        % type : string, 'box', 'slew', 'accum'. These types correspond to
        %        magnitude constraints, rate of change constraints, and
        %        integrated (accumulated) constraint, respectively. The
        %        'accum' option is useful to provide an actual magnitude
        %        constraint when the optimization is over \Delta u, rather
        %        than u.
        %        
        %        N.B. The 'accum' option will currently only work properly
        %        for a "one-shot" optimization, like time-optimal or CLQR.
        % 
        % bnds : is a 1x2 vector of upper and lower bounds. If type =
        % 'slew', then only a 1x1 double is required.
        
            if strcmp('box', type) & length(bnds)==1
                bnds = [-bnds, bnds];
            end

            if strcmp(type, 'box')
                self.ub = ones(self.N_mpc, 1)*bnds(2);
                self.lb = ones(self.N_mpc, 1)*bnds(1);
            elseif strcmp(type, 'slew')
                S = derMat(self.N_mpc);
                self.Ainq = [self.Ainq;
                             S];
                self.lbAinq_funs{end+1} = @(x0)repmat(-bnds(1), self.N_mpc-1, 1);
                self.ubAinq_funs{end+1} = @(x0)repmat(bnds(1), self.N_mpc-1, 1);
            elseif strcmp(type, 'accum')
               S = CondenCon.accumMat(self.N_mpc);
               self.Ainq = [self.Ainq; S];
               self.ubAinq_funs{end+1} = @(x0)repmat(bnds(2), self.N_mpc, 1);
               self.lbAinq_funs{end+1} = @(x0)repmat(bnds(1), self.N_mpc, 1);
                
            end
            
        end
        
    end % methods() ends here
    
    methods(Static)
        function A_mat = init_cond_resp_matrix(a, i_start, i_end, c)
        % Returns the matrix 
        %         | C*A^i_start   | 
        %         | C*A^i_start+1 | 
        % A_mat = | C*A^3         |
        %         |  :          |
        %         | C*A^i_end     |
        % 
        % such that if i_start = 0, i_end = N, then A_mat*x0 returns
        % the output trajectory (for zero input). If c is empty, then
        % c is set to the indentity, so that A_mat*x0 is the state
        % sequence. 
            if ~exist('c', 'var') || isempty(c)
                c = eye(size(a, 1));
            end
            A_mat = [];
            for i = i_start:i_end
                A_mat = [A_mat; c*a^i];
            end
            
        end
        
        function B_mat = zero_state_output_resp(sys, N, i_start)
        % Returns a matrix such that B_mat * [u0, ... u_N] yeilds the
        % forced zero state response. Nominally, 
        %
        %  B_mat = |0               |  <-- i_start = 0
        %          |CB 0            |  <-- i_start = 1
        %          |CAB CB 0        |
        %          |CA^2B CAB CB, 0 |
        % is the output of initCondRespMatrix.
        if ~exist('i_start', 'var')
            i_start = 0;
        end
            nu = size(sys.b, 2);
            no = size(sys.c, 1);
            ns = size(sys.b, 1);
            I = eye(ns);
            
            CB_col1 = [sys.c*I*0;
                       CondenCon.init_cond_resp_matrix(sys.a, 0, N-1, sys.c)]*sys.b;
            CB_col1 = CB_col1(i_start*no +1:end, :);
            
            block_size = [no, nu];
            B_mat = CondenCon.reverse_block_hankel(CB_col1, block_size);
        end
        
        function B_mat = zero_state_resp(sys, N, i_start)
        % Returns a matrix such that B_mat * [u0, ... u_N] yeilds the
        % forced zero state response. Nominally, 
        %
        %  B_mat = |  0          |  <-- i_start = 0
        %          |  B   0      |  <-- i_start = 1
        %          | AB   B  0   |
        %          |A^2B AB  B 0 |
        % is the output of initCondRespMatrix.
        if ~exist('i_start', 'var')
            i_start = 0;
        end
            nu = size(sys.b, 2);
            no = size(sys.c, 1);
            ns = size(sys.b, 1);
            I = eye(ns);
            
            B_col1 = [I*0;
                     CondenCon.init_cond_resp_matrix(sys.a, 0, N-1, I)]*sys.b;
            B_col1 = B_col1(i_start*ns +1:end, :);
            block_size = [ns, nu];
            B_mat = CondenCon.reverse_block_hankel(B_col1, block_size);
        end
        
        function H_mat = reverse_block_hankel(B_col1, block_size)
        % Given a matrix M partitioned as
        % 
        %      | M1 | 
        %      | M2 | 
        % M  = | M3 |
        %      |  : |
        %      | MN |
        %
        % where each block has size block_size(1) rows and block_size(2)
        % columns, the function returns
        %
        %  H_mat = |M1, 0,...     |
        %          |M2, M1, 0,    |
        %          |M3  M2, M1, 0 |
        %          |MN  M3, M2, M1|
        %
        % Not sure if there is a better name for this matrix but it is
        % somewhat similar to a hankel matrix, but with some things
        % reversed. 
            no = block_size(1);
            nu = block_size(2);
            H_mat = [];
            N = size(B_col1, 1)/no;
            for i=0:N-1 % Build accross columns
               H_mat = [H_mat, [zeros(no*i,nu); 
                               B_col1(1:end-(i)*no, :)] ];
            end
        end
        
        function M = accumMat(N)
            M = zeros(N, N);
            for k=1:N
               M(k, 1:k) = ones(1, k); 
            end
        end        
    end  % methods(Static) ends here.  
    
end

