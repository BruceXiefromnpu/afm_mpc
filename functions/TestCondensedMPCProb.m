classdef TestCondensedMPCProb < matlab.unittest.TestCase
% This class provides tests for the class CondenCon.
    properties
       sys;
       N_mpc;
       Q;
       R;
       S;
    end
    methods
        function self = TestCondensedMPCProb()
            % Data from the power amp system.
            A = [1.429,   0.7208;
                -0.4429,   0.2865];
            B = [-0.2916;
                0.779];
            C = [0.3434,  0.1286];
            Ts = 40e-6;
            self.sys = ss(A, B, C, 0, Ts);
            
            
            % Generate a positive definitate matrix.
            rng(1);
            Q = cov(randn(2,1), randn(2,1));
            self.Q = Q + eye(2);
            self.R = 10;
            self.S = randn(2,1);
        end
    end
    
    methods (Test)
        function test_build_mpc_problem(self)
            N = 4;
            A = self.sys.A;
            B = self.sys.B;
            R = self.R;
            Q = self.Q;
            S = self.S;
            S0 = S*0;
            
            I = eye(size(A,1));
            Z = B*0;
            AA = [A^0;
                 A^1;
                 A^2;
                 A^3;
                 A^4];
             
            BB = [Z,      Z     Z,   Z;
                  B,      Z,    Z,   Z;
                  A*B,    B,    Z,   Z;
                  A^2*B, A*B,   B,   Z;
                  A^3*B, A^2*B, A*B, B];
            
            RR = blkdiag(R, R, R, R);
            Qp = dare(A, B, Q, R, S0);
            QQ = blkdiag(Q, Q, Q, Q, Qp); %#ok<*PROP>

            H_expected = 2*(RR + BB'*QQ*BB);
            M_expected = 2*(AA'*QQ*BB)';
            
            [H, M] = CondensedMPCProb.build_mpc_problem(self.sys,...
                N, Q, R, Qp, S0);
            
            self.verifyEqual(H, H_expected, 'AbsTol', 1e-14);
            self.verifyEqual(M, M_expected, 'AbsTol', 1e-14);

            % Test with non-zero S
            Qp = dare(A, B, Q, R, S);
            QQ = blkdiag(Q, Q, Q, Q, Qp); %#ok<*PROP>
            SS = [blkdiag(S,S,S,S);
                  zeros(size(B,1), size(B,2)*N)];
            
            H_expected = 2*(RR + BB'*QQ*BB + BB'*SS + SS'*BB);
            M_expected = 2*(AA'*QQ*BB + AA'*SS)';
            [H, M] = CondensedMPCProb.build_mpc_problem(self.sys,...
                N, Q, R, Qp, S);
            
            self.verifyEqual(H, H_expected, 'AbsTol', 1e-14);
            self.verifyEqual(M, M_expected, 'AbsTol', 1e-14);
            
            % If unconstrained, we should get the same control as lqr, by
            % selection of Qp.
            K_lqr = dlqr(A, B, Q, R, S);
            K_equiv = [1, zeros(1, N-1)]* inv(H)*M;
            
            self.verifyEqual(K_equiv, K_lqr, 'AbsTol', 1e-14);

        end
        
        function test_accum_mat(self)
            rng(1);
            N = 10;
            x = randn(N);
            % There should be a better way to do this, but the only way to
            % test the non-static methods seems to be to create an instance
            % of a concrete class??
            mpcProb = condensedMPCprob_OA(self.sys, N, self.Q, self.Q, self.R);
            accM = mpcProb.accumMat(N);
            
            xsum_exp = cumsum(x);
            xsum = accM*x;
            self.verifyEqual(xsum, xsum_exp, 'AbsTol', 1e-15);
        end

        function test_discrete_lsim(self)
            N = 10;
            mpcProb = condensedMPCprob_OA(self.sys, N, self.Q, self.Q, self.R);
            
            x0 = [1;2];
            u = randn(N, 1);
            t = [0:1:N-1]'*self.sys.Ts;
            [~,~, x_exp] = lsim(self.sys, u, t, x0);
            x = mpcProb.discrete_lsim(u, x0);

            x_exp = [x_exp; (self.sys.A*x_exp(end,:)' + self.sys.B*u(end))'];
            
            self.verifyEqual(x, x_exp', 'AbsTol', 1e-15);
            
        end
 
        
    end



end

%%con.add_input_con('box')