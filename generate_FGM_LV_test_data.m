% clc
clear
%  close all

addpath('../functions')
% PATH_sim_model       = pwd;  % for simulink simulations

% ---- Paths for shuffling data to labview and back. ------
% ---------- Load Parametric Models  -----------

load(fullfile(PATHS.sysid, 'FRF_data_current_stage2.mat'))
whos
umax = 5;

TOL = .01;

%%
[Gvib, gdrift] = eject_gdrift(modelFit.models.G_uz2stage);
PLANT = ss(Gvib);
SYS = ss(Gvib);

SYS = balreal(SYS);
PLANT = balreal(SYS);
Nx = SSTools.getNxNu(PLANT);
T = diag(1./Nx)/10;
SYS = ss2ss(SYS, T);
PLANT = ss2ss(PLANT, T);

Nd = 9;
SYS.iodelay = 0;
SYS.InputDelay = 0;

PLANT.InputDelay = Nd;
PLANT = absorbDelay(PLANT);

SYS.InputDelay = Nd;
SYS = absorbDelay(SYS);
Ts  = SYS.Ts;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                  Design reference "trajectory"                          %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% We will track one setpoints. Number of samples for each setpoint is 800.
N1    = 800;
ref_f_1 = 6; % 1.5 to hit slew rate, 1.4 doesn't
t1 = [0:1:N1]'*Ts;
trun = Ts*(N1);
t = t1;
yref = [0*t1 + ref_f_1];
ref_0 = 0;
ref_traj = timeseries(yref, t);

%--------------------------------------------------------------------------
% 
% 3). Reduced order system for simulation.
sys_obs = absorbDelay(SYS);
Ns  = length(sys_obs.b);
% 3). LQR generation gain.
% -------------------------------------------------------------------------
% -------------------- Constrained LQR Stuff -----------------------------
du_max   = StageParams.du_max;

% Pull out open-loop pole-zero information.
[wp_real_x, wz_real_x] = w_zp_real(SYS);
rho_1 = wz_real_x(1)/wp_real_x(1);
rhos_x = [rho_1, 1., 1];

zeta_x = [.8, .7, .4, .4 .4];
gams_x = [1., 1., 1., 1., 1];

rad = 0.25;
pint = 0.8;
sys_recyc=SSTools.deltaUkSys(SYS);

P_x  = getCharDes(sys_recyc, gams_x, pint, zeta_x, rhos_x, rad);
[Chat, Dhat] = place_zeros(sys_recyc, P_x);
gam = 1;
Q1 = Chat'*Chat;
S1 = Chat'*Dhat;
R1 = Dhat'*Dhat + gam; % Plus gamma.

N_mpc = 12;
Qp = dare(sys_recyc.a, sys_recyc.b, Q1, R1, S1);
nu = 1;
mpcProb = condensedMPCprob_OA(sys_recyc, N_mpc, Q1, Qp, R1, S1);
Hmpc = mpcProb.H; Mmpc = mpcProb.M;
Ainq = [eye(nu*N_mpc);
        -eye(nu*N_mpc)];
binq = ones(2*N_mpc*nu, 1)*du_max;

% =========================================================================
% FGM in Floating Point
maxIter = 20;
fgm_fp = FGMprob_1(sys_recyc, N_mpc, Q1, Qp, R1, S1, du_max, maxIter);
I_H_mpc = fgm_fp.I_HL;
ML_x0  = fgm_fp.ML;
beta    = fgm_fp.beta;

z0 = ones(N_mpc*nu,1)*du_max;
y0 = ones(N_mpc*nu,1)*du_max;


fprintf('condition number of H: %f\n', mpcProb.kappa)
% ------------------------------------------------------------------------%
%                           Observer Gain                                 %
% ------------------------------------------------------------------------%
sys_obs = absorbDelay(SYS);
p_int_d = .85;
  % state disturbance does not work unless we put the deltaUk state in the
  % observer too. It probably could be made to, but I havnt worked that out.
Qw = sys_obs.b*sys_obs.b'*150;
Lx = dlqr(sys_obs.a', sys_obs.c', Qw, 1)';
[L_dist, sys_obsDist, IDENT_obs, eNs_12] = DistEst.output_dist_est(sys_obs,...
                                             Lx, p_int_d);
[Nx_r, Nx_d, Nu_r, Nu_d] = DistEst.steady_state_gains(sys_obs, sys_obs.b*0, 1);
[Nx] = SSTools.getNxNu(sys_recyc);
        
% ------------------------------------------------------------------------%
%                        Create FXP Representation
% ------------------------------------------------------------------------%
fgm_fxp.x_nw = 32;
fgm_fxp.x_nf = 27;

mpcp_use = fgm_fxp;
sim_struct = struct('mpcProb', mpcp_use, 'nw', nw, 'nf',nf, 'PLANT', PLANT);

nw = 32;
nf = 26;
nw_fgm = 32;
nf_fgm = 28;
fgm_fp = FGMprob_1(sys_recyc, N_mpc, Q1, Qp, R1, S1, du_max, maxIter);
fgm_fxp = FGMprob_fxp_1(sys_recyc, N_mpc, Q1, Qp, R1, S1, du_max,...
           maxIter, nw_fgm, nf_fgm);
         
A_obs_cl = sys_obsDist.a - L_dist*sys_obsDist.c;
fprintf('A_cl needs n_int = %d\n', ceil(log2(max(max(abs(A_obs_cl))))) + 1)
fprintf('L needs n_int = %d\n', ceil(log2(max(abs(L_dist)))) + 1)
fprintf('Nx needs n_int = %d\n', ceil(log2(max(abs(Nx)))) + 1)
fprintf('K needs n_int = %d\n', ceil(log2(max(abs(K_lqr)))) + 1)
fprintf('B needs n_int = %d\n', ceil(log2(max(abs(sys_obsDist.b))))+2)


sim_struct.du_max = fi(0.198, 1, 32, 27);
sim_struct.Nx_recyc = fi(Nx, 1, 32, 30);
sim_struct.L = fi(L_dist, 1, 32, 30);


sim_struct.sys_obs.a = fi(sys_obsDist.a -L_dist*sys_obsDist.c, 1, nw, nw-7);
sim_struct.sys_obs.b = fi(sys_obsDist.b, 1, nw, 29);
sim_struct.sys_obs.c = fi(sys_obsDist.c, 1, nw, 28);
sim_struct.x0_obs = sys_obsDist.b*0;

% MPP2 = condensedMPCprob_OA(sys_recyc, N_mpc, Q1, Qp, R1, S1);         
% MPP2.lb = zeros(N_mpc,1) - du_max;
% MPP2.ub = zeros(N_mpc,1) + du_max;

sim_struct.mpcProb = mpcp_use;
sim_struct.nw = nw;
sim_struct.nf = nf;


[y_mpc, U_mpc, U_nom, dU_mpc, Xhat_fxp, Xerr_fxp] = sim_AFM(sim_struct, ref_traj, 1);

mpcOpts = stepExpOpts('pstyle', '-r', 'TOL', TOL, 'y_ref', ref_f_1,...
                      'controller', mpcp_use, 'name',  'Simulation');

sim_exp = stepExpDu(y_mpc, U_mpc, dU_mpc, mpcOpts);
F1 = figure(58);
plot(sim_exp, F1, 'umode', 'both');
%%
F1 = figure(59); clf
subplot(3,1,1)
hold on, grid on;
plot(ref_traj.time, ref_traj.Data, '--k', 'LineWidth', .05);

%%
% =========================================================================
% FINALLY
% Save some data for labview testing.
kk = 5;
zi_dat = fgm_fxp.zi_mat(:, 1:kk);
yi_dat = fgm_fxp.yi_mat(:, 1:kk);
xerr_dat = Xerr_fxp.Data(1:kk, :)';
x_zi_yi_dat = [double(xerr_dat); double(zi_dat); double(yi_dat)];

MPCMatrix = [I_H_mpc, zeros(12,2), ML_x0, zeros(12, 2)];
mpcDataPath = fullfile(PATHS.labview, ['labview\UnitTests\fpga_harnesses',
                                  '\MPC_FXP_FGM_TestCase\fgm_matrix.csv']);

csvwrite(controlDataPath, MPCMatrix);
AllMatrix = packMatrixDistEstMPC_singleAxis(sys_obsDist, L_dist, Nx);

dpath = 'C:\Users\arnold\Documents\labview\qpPaper\UnitTests\fgm_result_data.csv';

fid = fopen(dpath, 'w+');
% fid = 1
fprintf(fid, '%f, %f, %f, %f, %f, %f\n', N_mpc, size(sys_recyc.b, 1),...
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









