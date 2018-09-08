clc
clear
%  close all

addpath('functions')
% PATH_sim_model       = pwd;  % for simulink simulations

% ---- Paths for shuffling data to labview and back. ------
% ---------- Load Parametric Models  -----------


umax = 5;

TOL = .01;

%%

[plants, frf_data] = CanonPlants.plants_ns14;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                  Design reference "trajectory"                          %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% We will track one setpoints. Number of samples for each setpoint is 800.
N1    = 400;
r1 = 2;
step_ref = StepRef([r1], N1);
yref = step_ref.yref;


% 3). LQR generation gain.
% -------------------------------------------------------------------------
% -------------------- Constrained LQR Stuff -----------------------------
du_max   = StageParams.du_max;

cmplx_rad = 0.9;
gam_mpc = 5;

[Q1, R0, S1, P_x] = build_control_constsigma(plants.sys_recyc, cmplx_rad);
  
N_mpc = 22;
Qp = dare(plants.sys_recyc.a, plants.sys_recyc.b, Q1, R0+gam_mpc, S1);

% ------------------------------------------------------------------------%
%                           Observer Gain                                 %
% ------------------------------------------------------------------------%
% sys_obs = absorbDelay(SYS);
% p_int_d = .85;
%   % state disturbance does not work unless we put the deltaUk state in the
%   % observer too. It probably could be made to, but I havnt worked that out.
% Qw = sys_obs.b*sys_obs.b'*150;
% Lx = dlqr(sys_obs.a', sys_obs.c', Qw, 1)';
% [L_dist, sys_obsDist, IDENT_obs, eNs_12] = DistEst.output_dist_est(sys_obs,...
%                                              Lx, p_int_d);

% ------------------------- Observer Gain ---------------------------------
can_obs_params = CanonObsParams_01();
can_obs_params.beta = 50;
[sys_obsDist, L_dist] = build_obs(plants.SYS, can_obs_params);

% 2). Design FeedForward gains.
[Nx, Nu] = SSTools.getNxNu(plants.sys_recyc);

        
% ------------------------------------------------------------------------%
%                        Create FXP Representation
% ------------------------------------------------------------------------%
fgm_fxp.x_nw = 32;
fgm_fxp.x_nf = 27;
maxIter = 20;

nw = 32;
nf = 26;
nw_fgm = 32;
nf_fgm = 28;
fgm_fp = FGMprob_1(plants.sys_recyc, N_mpc, Q1, Qp, R0+gam_mpc, S1, du_max, maxIter);
I_H_mpc = fgm_fp.I_HL;
ML_x0  = fgm_fp.ML;
beta    = fgm_fp.beta;

fgm_fxp = FGMprob_fxp_1(plants.sys_recyc, N_mpc, Q1, Qp, R0+gam_mpc, S1, du_max,...
           maxIter, nw_fgm, nf_fgm);
% fgm_fxp.ML = fi(fgm_fp.ML, 1, 32, 26);
fgm_fxp.x_nw = 32;
fgm_fxp.x_nf = 27;
         
A_obs_cl = sys_obsDist.a - L_dist*sys_obsDist.c;
fprintf('A_cl needs n_int = %d\n', ceil(log2(max(max(abs(A_obs_cl))))) + 1)
fprintf('L needs n_int = %d\n', ceil(log2(max(abs(L_dist)))) + 1)
fprintf('Nx needs n_int = %d\n', ceil(log2(max(abs(Nx)))) + 1)
fprintf('B needs n_int = %d\n', ceil(log2(max(abs(sys_obsDist.b))))+2)


du_max_fxp = fi(0.198, 1, 32, 27);
Nx_fxp = fi(Nx, 1, 32, 30);
L_fxp = fi(L_dist, 1, 32, 30);

sys_obs_fxp.a = fi(sys_obsDist.a -L_dist*sys_obsDist.c, 1, nw, nw-7);
sys_obs_fxp.b = fi(sys_obsDist.b, 1, nw, 29);
sys_obs_fxp.c = fi(sys_obsDist.c, 1, nw, 28);
x0_obs = sys_obsDist.b*0;

 sims_fxpm = SimAFM(plants.PLANT, fgm_fxp, Nx_fxp, sys_obs_fxp, L_fxp, du_max_fxp,...
    true, 'nw', nw, 'nf', nf);
[y_fxpm, U_full_fxpm, U_nom_fxpm, dU_fxpm, Xhat_fxpm, Xerr_fxpm] = sims_fxpm.sim(yref);
name = 'FXP MPC Sim. (unit-test)';
%
fxpm_Opts = stepExpDuOpts('pstyle', '--g', 'TOL', TOL, 'step_ref', step_ref,...
  'controller', fgm_fxp, 'name',  name);

sim_exp_fxpm = stepExpDu(y_fxpm, U_full_fxpm, dU_fxpm, fxpm_Opts);

F1 = figure(58);
plot(sim_exp_fxpm, F1, 'umode', 'both');
ax=subplot(3,1,1);
hold on, grid on;
plot(step_ref, ax,'--k', 'LineWidth', .05);


% =========================================================================
% FINALLY
% Save some data for labview testing.

fgm_fxp.reset_warm_start_data();

for k = 1:10
  fgm_fxp.solve(Xerr_fxpm.Data(k,:)');
end


%%
kk = 5;
zi_dat = fgm_fxp.zi_mat(:, 1:kk);
yi_dat = fgm_fxp.yi_mat(:, 1:kk);
xerr_dat = Xerr_fxpm.Data(1:kk, :)';
x_zi_yi_dat = [double(xerr_dat); double(zi_dat); double(yi_dat)];

MPCMatrix = [I_H_mpc, zeros(N_mpc,2), ML_x0, zeros(N_mpc, 2)];
fname = sprintf('fgm_matrixN%d.csv', N_mpc)

lv_test_root = fullfile(PATHS.labview, 'UnitTests', 'fpga_harnesses',...
                                  'MPC_FXP_FGM_TestCase')
mpcDataPath = fullfile(lv_test_root, fname);
csvwrite(mpcDataPath, MPCMatrix);


fname_results = sprintf('fgm_result_data_N%d.csv', N_mpc)
dpath = fullfile(lv_test_root, fname_results);

fid = fopen(dpath, 'w+');
% fid = 1
fprintf(fid, '%f, %f, %f, %f, %f, %f\n', N_mpc, size(plants.sys_recyc.b, 1),...
  -beta, beta+1, du_max, maxIter);
for k=1:size(x_zi_yi_dat, 1)
  fprintf(fid, '%.12f,', x_zi_yi_dat(k,:));
  fprintf(fid, '\n');
end

fclose(fid);

% 
% M_fxp = fgm_fxp.ML;
% xx = double(xerr_dat(:,1));
% xx = fi(xx, 1, fgm_fxp.x_nw, fgm_fxp.x_nf);
% f = fi(M_fxp*xx, 1, fgm_fxp.nw, fgm_fxp.nf);
% idx = 2;
% z0 = zi_dat(:,idx);
% y0 = zi_dat(:,idx);
% fgm_fxp.warm_start_data = fi(zeros(24, 1), 1, 32, 28);
% 
% [z1, y1] = fgm_fxp.call_qp_solver(f, z0, y0)
% 
% [zi_dat(:,2), z1]
% 









