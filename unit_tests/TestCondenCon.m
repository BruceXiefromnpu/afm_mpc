classdef TestCondenCon < matlab.unittest.TestCase
% This class provides tests for the class CondenCon.
    properties
       sys;
       N_mpc;
    end
    methods
        function self = TestCondenCon()
            % Data from the power amp system.
            A = [1.429,   0.7208;
                -0.4429,   0.2865];
            B = [-0.2916;
                0.779];
            C = [0.3434,  0.1286];
            Ts = 40e-6;
            self.sys = ss(A, B, C, 0, Ts);
        end
    end
    
    methods (Test)
        function test_add_input_con(self)
        % Verify that input constraints are handled properly.
            con = CondenCon(self.sys, [0;0], 4);
            
            % These should get added to lb and ub
            con.add_input_con('box', [-0.4, 0.4])
            % These should get added to function handles.
            con.add_input_con('slew', [0.5])
            self.verifyEqual(length(con.ubAinq_funs), 1)
            self.verifyEqual(length(con.lbAinq_funs), 1)
            
            con.update_binq()
            lbAinq_exp = [-0.5; -0.5; -0.5];
            ubAinq_exp = [0.5; 0.5; 0.5];
            lb_exp = [-0.4; -0.4; -0.4; -0.4];
            ub_exp = [0.4; 0.4; 0.4; 0.4];
            Ainq_exp = [-1, 1, 0, 0;
                         0, -1, 1, 0;
                         0,  0, -1, 1];
            self.verifyEqual(con.Ainq, Ainq_exp)
            
            self.verifyEqual(con.lbAinq, lbAinq_exp)
            self.verifyEqual(con.ubAinq, ubAinq_exp)
            
            self.verifyEqual(con.lb, lb_exp)
            self.verifyEqual(con.ub, ub_exp)

        end
        function test_update_sys(self)
        % Verify that we can propograte the system properly. 
            x0 = [0;0];
            uk = 5.3;
            con = CondenCon(self.sys, x0, 4);
            con.update_sys(uk);
            % A*x0 = 0, when x0 = 0, so xk+1 = b*uk
            self.verifyEqual(con.x0, self.sys.b*uk)
        end
        
        function test_add_state_con(self)
        % Verify that state constraints get added properly.
            x0 = [0;0];
            N = 4;
            yMax = 0.6;

            con = CondenCon(self.sys, x0, N);
            con.add_state_con('box', yMax);
            con.update_binq()
            % Pull out for easier access.
            A = self.sys.a;
            B = self.sys.b;
            C = self.sys.c;
            I = eye(size(B, 1));
            f = [C;
                 C*A^1;
                 C*A^2;
                 C*A^3;
                 C*A^4];
            Ainq_exp = [0,   0  0, 0;
                        C*B, 0, 0, 0;
                        C*A*B, C*B, 0, 0;
                        C*A^2*B, C*A*B, C*B, 0;
                        C*A^3*B, C*A^2*B, C*A*B, C*B];
            
            self.verifyEqual(con.Ainq, Ainq_exp, 'AbsTol', 1e-15);
            % We start with x0=0, so expect that lb(ub)Ainq = \pm 0.6
            con.update_binq();
            ubAinq_exp = ones(N+1, 1)*yMax;
            self.verifyEqual(con.ubAinq, ubAinq_exp);
            self.verifyEqual(con.lbAinq, -ubAinq_exp);
            
            % Now, make x0 non zero, and test again.
            con.x0 = [1;2];
            con.update_binq();
            ubAinq_exp = ones(N+1, 1)*yMax - f*con.x0;
            lbAinq_exp = -ones(N+1, 1)*yMax - f*con.x0;
            self.verifyEqual(con.ubAinq, ubAinq_exp, 'AbsTol', 1e-15);
            self.verifyEqual(con.lbAinq, lbAinq_exp, 'AbsTol', 1e-15);
            
        end
        
        function test_add_state_and_input_con(self)
            con = CondenCon(self.sys, [0;0], 4);

            con.add_state_con('box', 0.6);

            con.add_input_con('box', [-0.4, 0.4])
            con.add_input_con('slew', [0.4])
            con.update_binq()

        end
        function test_initCondRespMatrix(self)
        % Test the initial condition response matrix generation 
           A = [1, 2;
               3, 4];
           C = [1, 2];
           
           A_initcond_expected = [C;
                                  C*A;
                                  C*A^2];
          A_initcond = CondensedTools.init_cond_resp_matrix(A, 0, 2, C);
          
          self.verifyEqual(A_initcond_expected, A_initcond, 'AbsTol', 1e-15)
            
        end
    
        function test_reverse_block_hankel(self)
        % Verify that we construct the "reverse" block hankel matrix as
        % expected. 
           A = [1, 2;
               3, 4];
           B = [1;2];
           I = eye(2);
           Z = I*B*0;

          Bmat_col1 = [Z;
                       B;
                       A*B];
          Bmat_exp = [Z, Z,   Z;
                      B,  Z,  Z;
                      A*B, B, Z];
          B_mat = CondensedTools.reverse_block_hankel(Bmat_col1, [2, 1]);

          self.verifyEqual(Bmat_exp, B_mat, 'AbsTol', 1e-15)

        end
        
        function test_zero_state_resp(self)
           A = [1, 2; 
               3, 4];
           B = [1;2];
           C = [1 0];
           G = ss(A, B, C, 0, 0.01);
           I = eye(2);
           Z = I*B*0;

          Bmat_exp = [Z,  Z,  Z;
                      B,  Z,  Z;
                      A*B, B, Z];
          B_mat = CondensedTools.zero_state_resp(G, 2);

          self.verifyEqual(Bmat_exp, B_mat, 'AbsTol', 1e-15)

        end
        
        function test_zero_state_output_resp(self)
           A = [1, 2; 3, 4];
           B = [1;2];
           C = [1 0];
           G = ss(A, B, C, 0, 0.01);
           I = eye(2);
           Z = C*I*B*0;

          Bmat_exp = [Z,  Z, Z;
                      C*B,  Z, Z;
                      C*A*B, C*B, Z];
          B_mat = CondensedTools.zero_state_output_resp(G, 2);

          self.verifyEqual(Bmat_exp, B_mat, 'AbsTol', 1e-15)
          
          % Check it works with i_start
          Bmat_exp = [C*B,  Z;
                     C*A*B, C*B];
          B_mat = CondensedTools.zero_state_output_resp(G, 2, 1);

          self.verifyEqual(Bmat_exp, B_mat, 'AbsTol', 1e-15)

        end        
        
    end



end

%%con.add_input_con('box')