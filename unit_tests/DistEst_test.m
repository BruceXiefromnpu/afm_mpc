
function tests = DistEst_test()
  tests = functiontests(localfunctions);
end


function testDistEst_state_construction(tc)
  % Check that the physical construction is correct
  A = [1, 2;
       0  1];
  B = [1; 1];
  C = [1, 1];
  
  sys = ss(A, B, C, 0, 0.1);
  [~, sys_] = DistEst.state_dist_est(sys, [1;1], 1, B, 1);
  
  A_ = [A, B; 0, 0, 1];
  B_ = [B; 0];
  C_ = [C, 0];
  %keyboard
  tc.verifyEqual(sys_.a, A_);
  tc.verifyEqual(sys_.b, B_);
  tc.verifyEqual(sys_.c, C_);
  
  % Check it works if we dont provide Ad and Bd
  [~, sys_] = DistEst.state_dist_est(sys, [1;1], 1);
  
  % expect Bd = [0;1]
  A_ = [A, [0;1];
        0, 0, 1];
  B_ = [B; 0];
  C_ = [C, 0];
  %keyboard
  tc.verifyEqual(sys_.a, A_);
  tc.verifyEqual(sys_.b, B_);
  tc.verifyEqual(sys_.c, C_);  
  
end



function testDistEst_state(testCase)
load(fullfile(PATHS.sysid, 'FRF_data_current_stage2.mat'));

sys = ss(modelFit.models.G_uz2stage);
Q1 = sys.b*sys.b'*200;

Lx = dlqr(sys.a', sys.c', Q1, 1)';

sys_cl = ss(sys.a - Lx*sys.c, sys.b, sys.c, 0, sys.Ts);


p_int = 0.8;
[L, sys_dist_obs, Iobs, ens1] = DistEst.state_dist_est(sys, Lx, p_int);

sys_cl_2 = ss(sys_dist_obs.a - L*sys_dist_obs.c, sys_dist_obs.b, sys_dist_obs.c, 0, sys.Ts);

p1 = sort([pole(sys_cl); p_int]);
p2 = sort(pole(sys_cl_2));

TOL = 1e-10;
testCase.verifyEqual(p1, p2, 'AbsTol', TOL);

% Check the selector matrices work
% test vector
rng(1);
x = randn(size(sys.b,1)+1, 1);

testCase.verifyEqual(x(1:end-1), Iobs*x, 'AbsTol', 1e-15)
testCase.verifyEqual(x(end), ens1*x, 'AbsTol', 1e-15)



end

function testDistEst_state_deltaUk(testCase)
  % Check things work with a pole at z=1.
load(fullfile(PATHS.sysid, 'FRF_data_current_stage2.mat'));

sys = SSTools.deltaUkSys(ss(modelFit.models.G_uz2stage));
Q1 = sys.b*sys.b'*200;

Lx = dlqr(sys.a', sys.c', Q1, 1)';

sys_cl = ss(sys.a - Lx*sys.c, sys.b, sys.c, 0, sys.Ts);


p_int = 0.8;
[L, sys_dist_obs] = DistEst.state_dist_est(sys, Lx, p_int);

sys_cl_2 = ss(sys_dist_obs.a - L*sys_dist_obs.c, sys_dist_obs.b, sys_dist_obs.c, 0, sys.Ts);

p1 = sort([pole(sys_cl); p_int]);
p2 = sort(pole(sys_cl_2));

TOL = 1e-10;
testCase.verifyEqual(p1, p2, 'AbsTol', TOL);

end

function testDistEst_output_deltaUk(testCase)
  % Check things work with a pole at z=1.
load(fullfile(PATHS.sysid, 'FRF_data_current_stage2.mat'));

sys = ss(modelFit.models.G_uz2stage);
% SSTools.deltaUkSys(ss(modelFit.models.G_uz2stage));
Q1 = sys.b*sys.b'*200;

Lx = dlqr(sys.a', sys.c', Q1, 1)';

sys_cl = ss(sys.a - Lx*sys.c, sys.b, sys.c, 0, sys.Ts);


p_int = 0.8;
[L, sys_dist_obs] = DistEst.output_dist_est(sys, Lx, p_int);

sys_cl_2 = ss(sys_dist_obs.a - L*sys_dist_obs.c, sys_dist_obs.b, sys_dist_obs.c, 0, sys.Ts);

p_exp = sort([pole(sys_cl); p_int]);
p_act = sort(pole(sys_cl_2));

TOL = 1e-10;
testCase.verifyEqual(p_act, p_exp, 'AbsTol', TOL);

end

function testDistEst_steady_state_gains(tc)
  load(fullfile(PATHS.sysid, 'FRF_data_current_stage2.mat'));

  sys = ss(modelFit.models.G_uz2stage);
  Bd = sys.b;
  [Nx_r, Nx_d, Nu_r, Nu_d] = DistEst.steady_state_gains(sys, Bd);
  
  ns = size(sys.b,1);
  M = [eye(ns) - sys.a, -sys.b;
     sys.c,      0];
  dss = 1;
  rss = 2;
  
  xss = Nx_d*dss + Nx_r*rss;
  uss = Nu_d*dss + Nu_r*rss;
  
  lhs = M*[xss; uss];
  rhs = [Bd; 0]*dss + [Bd*0; 1]*rss;
  
  tc.verifyEqual(lhs, rhs, 'AbsTol', 1e-12);
end





