
addpath('~/matlab/afm_mpc_journal/functions')
addpath('~/matlab/afm_mpc_journal/functions/canon')
plants = CanonPlants.plants_with_drift_inv(true);

can_cntrl = CanonCntrlParams_01(plants.sys_recyc);
[Q1, R0, S1] = build_control(plants.sys_recyc, can_cntrl);

gam_lin = 3;
gam_mpc = 1;
R1 = R0 + gam_mpc;

K_lqr = dlqr(plants.sys_recyc.a, plants.sys_recyc.b, Q1, R0+gam_lin, S1);
sys_cl = SSTools.close_loop(plants.sys_recyc, K_lqr);

can_obs_params = CanonObsParams_01();
[sys_obsDist, L_dist] = build_obs(plants.SYS, can_obs_params);

f10 = figure(10); clf
pzplotCL(sys_cl, K_lqr, [], f10);


% -------------- Plot Observer Closed-loop poles -------------
F20 = figure(20); clf
p_plant = pole(plants.PLANT);
z_plant = zero(plants.PLANT);
hpobs = plot(real(p_plant), imag(p_plant), 'xk', 'MarkerSize', 8);
hpobs.DisplayName = 'System Pole';

hold on
hzobs = plot(real(z_plant), imag(z_plant), 'ok', 'MarkerSize', 8);
hzobs.DisplayName = 'System Zero';    

hold on
opts.pcolor = 'r';

[~, hpcl] = pzplotCL(sys_obsDist, [], L_dist, gcf, opts, 'MarkerSize', 8);
hpcl.DisplayName = 'C.L Obs Pole';

xlim([-0.05, 1.05])
ylim([-0.4, 0.4])
title('');
xlabel('Re')
ylabel('Im')
leg1 = legend([hpobs, hzobs, hpcl]);
set(leg1, 'location', 'SouthWest', 'box', 'off', 'FontSize', 14);

if saveon
  saveas(F20, fullfile(PATHS.jfig, 'obs_cl.svg'))
end

    
