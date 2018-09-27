classdef FGMprob_1 < CondensedMPCProb
% Holds problem data for the Fast gradient Method mpc QP solver.
% Properties:
%         I_HL;
%         ML;
%         beta;
%         beta_1;
%         L;
%         mu;
%         uMax;
%         condition;
%
% Construction:
% FGM(mpcProb, uMax)
%   mpcProb: an condensedMPCprob object
%   uMax: saturation bound. Optional.

    properties
        I_HL; % I - H/L, where H is the problem hessian and L is the max ev.
        ML; % M/L
        beta; % (sqrt(L) - sqrt(mu))/(sqrt(L) + sqrt(mu));
        beta_1; % 1-beta.
        L; % Max eigenvalue
        mu; % Min eigenvalue
        uMax; % upper/lower bound for the box constraint.
        maxIter; % Number of solver iterations for the fast gradient method.
        n_warmstart;
    end

    methods
        function self = FGMprob_1(sys, N, Q,Qp, R, S, uMax, maxIter)
          % fgm_prob = FGMprob_1(sys, N, Q,Qp, R, S, uMax, maxIter)
          self.N_mpc = N;
          self.nu = size(sys.b,2);
          [H, M] = CondensedMPCProb.build_mpc_problem(sys,N, Q, R, Qp, S);
          
          eigs = eig(H);
          L    = max(eigs);
          mu   = min(eigs);

          if mu < 0
            error('Hessian must be positive definite. Smallest eigenvalue is mu = %f', mu)
          end
          self.kappa = cond(H);
          self.H = H;
          self.M = M;
          self.I_HL   = eye(N) - H/L;
          self.ML     = M/L;
          self.beta   = (sqrt(L) - sqrt(mu))/(sqrt(L) + sqrt(mu));
          self.beta_1 = self.beta + 1;
          self.mu     = mu;
          self.L      = L;
          n_warmstart = 2*self.N_mpc*self.nu;
          self.n_warmstart = n_warmstart;
          %            obj.n_warmstart = 2*obj.N_mpc*obj.nu;
          self.warm_start_data = zeros(n_warmstart,1);
          %            obj.condition = mpcProb.condition;
          self.maxIter = maxIter;
          if exist('uMax', 'var')
            self.uMax = uMax;
          end

        end

        function reset_warm_start_data(self)
        % reset_warm_start_data(self, ziyi)
        % Resets the warm start data. If ziyi is provided, will
        % reset self.warm_start_data = ziyi. Otherwise, set it to
        % all zeros.
            self.warm_start_data = zeros(self.n_warmstart,1);
        end

        function [uk, ziyi] = solve(self, xk_1)
        % Solves via the fast gradient method a QP given an instance of
        % FGMprob, an initial condition xk_1, and warmstart/initial condition
        % data ziyi = [zi; yi]
        % NOTE: This method should be overloaded from the Parent
        % class, CondensedMPCProb.
              
            nu = self.nu;
            ziyi = self.warm_start_data;
            z_i = ziyi(1:self.n_warmstart/2);
            y_i = ziyi(self.n_warmstart/2+1:end);
            % shift zi and yi for warm start
            z_i = [z_i(nu+1:end); z_i(end-nu+1:end)];
            y_i = [y_i(nu+1:end); y_i(end-nu+1:end)];

            f = self.ML*xk_1;
            [z_i, y_i] = self.call_qp_solver(f, z_i, y_i);

            ziyi = [z_i; y_i];
            self.warm_start_data = ziyi;
            uk = z_i(1:self.nu);
        end
          
        function [z_i, y_i] = call_qp_solver(self, f, z0, y0)
        % [z_i, y_i] = call_qp_solver(obj, f, z0, y0)
        % mu          = min(eigs);
        % beta        = (sqrt(Lip) - sqrt(mu))/(sqrt(Lip) + sqrt(mu));
        % I_H3        = eye(N_mpc) - H3_mpc/Lip;
        % H2_L_x0_mpc = H2_x0_mpc/Lip;
        % H2_L = H2_L_x0_mpc*(x0-xss*uss)        
          N = length(f);
          % init:
          z_i   = z0;
          z_i_1 = z_i;
          y_i   = y0;
          du_max = self.uMax;
          for i = 1:self.maxIter
            % gradient evauluation:
            t_i = self.I_HL*y_i - f;
            % projection:
            z_i_1 =  max(min(t_i, du_max), -du_max);
            % update.
            y_i = (1+self.beta)*z_i_1   -self.beta*z_i;
            z_i = z_i_1;
          end
          z_i = z_i_1;
        end

    end % methods

end


% mu          = min(eigs);
% beta        = (sqrt(Lip) - sqrt(mu))/(sqrt(Lip) + sqrt(mu));
% I_H3        = eye(N_mpc) - H3_mpc/Lip;
% H2_L_x0_mpc = H2_x0_mpc/Lip;
% H2_L = H2_L_x0_mpc*(x0-xss*uss)






