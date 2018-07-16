% This script plots the LQR based root locus for the "lowgain"
% version of things.


saveon = true;
plants = CanonPlants.plants_ns14();
sys_recyc = plants.sys_recyc;
Ts = sys_recyc.Ts;

% can_cntrl = CanonCntrlParams_01(plants.SYS);
can_cntrl = CanonCntrlParams_ns14();
[Q1, R0, S1] = build_control(sys_recyc, can_cntrl);


f3 = figure(3); clf
P_x  = getCharDes(sys_recyc, can_cntrl.gam_s, can_cntrl.pint,...
                             can_cntrl.zeta_s, can_cntrl.rho_s, can_cntrl.rad);
plot(real(P_x), imag(P_x), 'ob')
hold on
[ax, C_hand] = lqr_locus(sys_recyc, Q1, 1, S1, .001, 1000);
C_hand.Label.Interpreter = 'latex';
C_hand.Label.FontSize = 14;
C_hand.Location = 'eastoutside';
C_hand.Position = [0.8406    0.1121    0.0473    0.8715];
ax.Position =  [0.0862 0.1123 0.7369 0.8713];
C_hand.Label.String = '$\gamma$';
C_hand.Label.Units = 'normalized';
C_hand.Label.Rotation = 0;
C_hand.Label.Position = [2.5214    0.5549         0];

xlim([-0.3, 1])
ylim([-0.35, 0.35])
% C_hand.Label.String = '$R_o + \gamma$';

xlabel('Re')
ylabel('Im')


if saveon
  saveas(f3, fullfile(PATHS.jfig, 'lqr_locus.svg'))
end
