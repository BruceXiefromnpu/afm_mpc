TOL = .01;
tol_mode = 'abs';
plotstate = false;
plotpoles = false;
plotdriftpoles = false;
plotdriftbode = false;
saveon = false;

md = 1;
% ------- Load Plants -----
with_hyst = true;
%%% plants = CanonPlants.plants_with_drift_inv(with_hyst);

% [plants, frf_data] = CanonPlants.plants_drift_inv_hyst_sat();
[plants, frf_data] = CanonPlants.plants_ns14;
% plants.PLANT = plants2.PLANT; %*dcgain(plants.Gvib)/dcgain(plants2.PLANT);
% plants.gdrift = plants2.gdrift;

Ts  = plants.SYS.Ts;
if md == 2
  plants.gdrift = zpk([], [], 1, Ts);
  plants.gdrift_inv = zpk([], [], 1, Ts);
end
plants.gdrift = plants.gdrift_1p0;
plants.gdrift_inv = 1/plants.gdrift_1p0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                  Design reference "trajectory"                          %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Get a ref trajectory to track.
N    = 8000;
r1 = 1;
trajstyle =1;
if trajstyle == 1
  step_ref = StepRef(r1, N);
  yref = step_ref.yref;
  dist_traj = yref;
  dist_traj.Data = dist_traj.Data*0 + 0.0;
  yref.Data = yref.Data*0;
end
rw = 2e-07;
rng(1);
thenoise = timeseries(mvnrnd(0, rw, length(yref.Time))*1, yref.Time);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
% Design control/estimator gains. This gains are what we actually         %
% use in both the simulation and experimentant system.                    %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 3). LQR generation gain.
% -------------------------------------------------------------------------
% -------------------- Constrained LQR Stuff ------------------------------

% Adjust the du_max to account for the gain of gdrift_inv.
du_max_orig = StageParams.du_max;
if md == 1
  du_max = du_max_orig/norm(plants.gdrift_inv, Inf);
else
  du_max = du_max_orig;
end
% du_max = 1000;

can_cntrl = CanonCntrlParams_ns14(plants.SYS);
% can_cntrl = can_cntrl.aggressive_params();
[Q1, R0, S1] = build_control(plants.sys_recyc, can_cntrl);
gam_lin = 30;
gam_mpc = .2;

R1 = R0 + gam_mpc;
% gams = linspace(0.0001, 200, 200);
gams = logspace(log10(0.1), log10(100), 50);
% gams = .1;
sixsigs = gams*0;
kappas = gams*0;

p_int = 0.8; cmplx_rad = 0.9; rho_s = [2, 1]; rad = 0.5;
Px = getCharDes_const_sig(plants.sys_recyc, p_int, cmplx_rad, rho_s, rad);
[Chat, Dhat] = place_zeros(plants.sys_recyc, Px);
Q1 = Chat'*Chat;
S1 = Chat'*Dhat;
R0 = Dhat'*Dhat;
 
% gam_lin = 12;
% gam_mpc = 2;
% R1 = R0 + gam_mpc;

N_mpc = 12;

for k=1:length(gams)

% feedback gain   
K_lqr = dlqr(plants.sys_recyc.a, plants.sys_recyc.b, Q1, R0+gams(k), S1);
Qp = dare(plants.sys_recyc.a, plants.sys_recyc.b, Q1, R0+gams(k), S1);

mpcProb = condensedMPCprob_OA(plants.sys_recyc, N_mpc, Q1, Qp, R0+gams(k), S1);
kappas(k) = mpcProb.kappa;

% ------------------------- Observer Gain ---------------------------------
can_obs_params = CanonObsParams_01();
can_obs_params.beta = 100;
[sys_obsDist, L_dist] = build_obs(plants.SYS, can_obs_params);


[Sens, Hyd, Hyr, Hyeta] = ss_loops_delta_dist(plants.SYS, plants.sys_recyc, sys_obsDist, K_lqr, L_dist);

Ry_est = cov(y_lin_fp_sim.Data);
Pcl = dlyap(Hyeta.a, Hyeta.b*rw*Hyeta.b');
Ry_calc = Hyeta.c*Pcl*Hyeta.c' + 0*rw;

sixsig = 6*sqrt(Ry_calc);
needed_14volt_cov = 1/512;
% fprintf('six-sigma cov: %.4g\n14v x 512 pix: %.4g\n', sixsig, needed_14volt_cov);

sixsigs(k) = sixsig;

end
%%
figure(21)
plot(gams, sixsigs)

xlabel('$\gamma$')
ylabel('6-$\sigma$')
%
sixsigpix_mu = sixsigs*512*5;

figure(22); clf
yyaxis left
plot(gams, sixsigpix_mu)
ylim([0, max(sixsigpix_mu)])
plotshaded(gams, [sixsigpix_mu; gams*0+max(sixsigpix_mu)], 'b', .2);
xlabel('$\gamma$')
ylabel('range [$\mu$m] (512 pixels)')
ax = gca();

% ax.YTick = [0:2.5:20]
grid on
yyaxis right
plot(gams, kappas)
ylabel('$\kappa$')






