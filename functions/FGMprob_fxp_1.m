classdef FGMprob_fxp_1 < CondensedMPCProb
% Holds problem data for the Fast gradient Method mpc QP solver.
% Properties:
%         I_HL;
%         ML;
%         beta;
%         beta_1;
%         L;
%         mu;
%         uMax;
%         kappa;
%         maxIter;
%         nw;
%         nf;
%         ziyi;
% Construction:
% FGM(mpcProb, uMax)
%   mpcProb: an condensedMPCprob object
%   uMax: saturation bound. Optional.

    properties
      I_HL;
      ML;
      %kappa;
      beta;
      beta_1;
      L;
      mu;
      uMax;
      maxIter;
      nw;
      nf;
      %%%%%%%%%%%%%%%%%
      n_warmstart;

      x_nw;
      x_nf;
      zi_mat;
      yi_mat;
    end

    methods
        function self = FGMprob_fxp_1(sys, N, Q,Qp, R, S, uMax, maxIter, nw, nf)

          self.N_mpc = N;
          self.nu = size(sys.b,2);
          self.ns = size(sys.b,1);
          [H, M] = CondensedMPCProb.build_mpc_problem(sys,N, Q, R, Qp, S);
          eigs = eig(H);
          L    = max(eigs);
          mu   = min(eigs);

           if mu < 0
              error(['Hessian must be positive definite. Smallest ',...
                     'eigenvalue is mu = %f'], mu)
           end
           self.kappa = cond(H);
           self.I_HL   = fi(eye(N) - H/L, 1, nw, nf);
           self.ML     = fi(M/L, 1, nw, nf);
           
           beta_fp = (sqrt(L) - sqrt(mu))/(sqrt(L) + sqrt(mu));
           self.beta   = fi(beta_fp, 1, nw, nf);
           self.beta_1 = self.beta + fi(1,1,2,0);
           self.mu     = mu;
           self.L      = L;
           self.n_warmstart = 2*self.N_mpc*self.nu;
           self.warm_start_data = [];
           self.maxIter = maxIter;

           if isa(uMax, 'embedded.fi')
             self.uMax = uMax;
           else
             self.uMax = fi(uMax, 1, nw, nf);
           end

           ziyi = zeros(self.n_warmstart,1);
           self.warm_start_data = fi(ziyi, 1, self.uMax.WordLength,...
                                     self.uMax.FractionLength);
           self.nw = nw;
           self.nf = nf;

        end

        function reset_warm_start_data(self, ziyi)
        % reset_warm_start_data(self, ziyi)
        % Resets the warm start data. If ziyi is provided, will
        % reset self.warm_start_data = ziyi. Otherwise, set it to
        % all zeros. 
          if exist('ziyi', 'var')
            self.warm_start_data = ziyi;
          else
            self.warm_start_data = fi(zeros(self.n_warmstart,1), 1, ...
                                      self.nw, self.nf);
          end
          n_ws = self.n_warmstart;
          self.zi_mat = self.warm_start_data(1:n_ws/2);
          self.yi_mat = self.warm_start_data(n_ws/2+1:end);;
        end

        function [uk] = solve(self, xk_1)
        % [uk] = solve(self, xk_1)
        % Solve the MPC problem (via FGM, obviously) given the
        % starting condition xk_1.
        % NOTE: This method should be overloaded from the Parent
        % class, CondensedMPCProb.
            %#codegen
            xk_1 = fi(xk_1, 1, self.x_nw, self.x_nf);
            nu = self.nu;
            z_i = self.warm_start_data(1:self.n_warmstart/2);
            y_i = self.warm_start_data(self.n_warmstart/2+1:end);
            % shift zi and yi for warm start

            z_i = [z_i(nu+1:end); z_i(end-nu+1:end)];
            y_i = [y_i(nu+1:end); y_i(end-nu+1:end)];
            %keyboard
            f = fi(self.ML*xk_1, 1, self.nw, self.nf);
            [z_i, y_i] = self.call_qp_solver(f, z_i, y_i);
            self.zi_mat = [self.zi_mat, z_i(:)];
            self.yi_mat = [self.yi_mat, y_i(:)];
            
            self.warm_start_data = [z_i; y_i];
            uk = double(z_i(1:self.nu));

        end
          
        function [z_i, y_i] = call_qp_solver(self, f, z0, y0)
        % [z_i, y_i] = call_qp_solver(obj, f, z0, y0)
        % mu          = min(eigs);
        % beta        = (sqrt(Lip) - sqrt(mu))/(sqrt(Lip) + sqrt(mu));
        % I_H3        = eye(N_mpc) - H3_mpc/Lip;
        % H2_L_x0_mpc = H2_x0_mpc/Lip;
        % H2_L = H2_L_x0_mpc*(x0-xss*uss)
          
          N = length(f);
          I_HL = self.I_HL;
          nw = self.nw;
          nf = self.nf;
          du_max = self.uMax;
          beta = self.beta;
          % init:
          z_i   = fi(z0, 1, nw, nf);
          z_i_1 = z_i;
          y_i   = fi(y0, 1, nw, nf);

          for i = 1:self.maxIter
            % gradient evauluation:
            t_i = fi(I_HL*y_i - f, 1, nw, nf);
            % Projection:
            z_i_1 =  max(min(t_i, du_max), -du_max);
            % update.
            y_i = accumpos(fi((1+beta)*z_i_1, 1, nw, nf),   fi(-beta*z_i, 1, nw, nf));
            z_i = z_i_1;
          end
          z_i = z_i_1;
        end % call_qp_solver
    end % methods

end



