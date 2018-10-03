% This script plots the LQR based root locus for the "lowgain"
% version of things.


saveon = true;
lqr_locus_figfile = 'lqr_locus_choozezet.svg';

plants = CanonPlants.plants_ns14(9, 2);
sys_recyc = plants.sys_recyc;
Ts = sys_recyc.Ts;


can_cntrl = CanonCntrlParams_ns14();
[Q1, R0, S1] = build_control(sys_recyc, can_cntrl);


f3 = figure(3); clf
P_x  = getCharDes(sys_recyc, can_cntrl.gam_s, can_cntrl.pint,...
                             can_cntrl.zeta_s, can_cntrl.rho_s, can_cntrl.rad);
plot(real(P_x), imag(P_x), 'ok')
hold on
ax = gca();

rgb1 = [0.230, 0.299, 0.754];
rgb2 = [0.706, 0.016, 0.150];
s_ = linspace(0,1, 400);
color_map = diverging_map(s_, rgb1, rgb2);

color_map = diverging_map(s_, rgb1, rgb2);

[ax] = lqr_locus(sys_recyc, Q1, 1, S1, .001, 1000, ax,...
  'color_map', color_map);

caxis([gam_min, gam_max]);
C_hand = colorbar('Location', 'eastoutside'); 
tck_labs = [0.001, 0.01, 0.1, 1, 10, 100, 1000, 10000];

C_hand.Ticks = linspace(gam_min, gam_max, length(tck_labs));
C_hand.TickLabels = tck_labs;


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
    fpath = fullfile(PATHS.jfig, lqr_locus_figfile);
  fprintf('saving Figure number %d as\n %s\n', f3.Number, fpath);
  saveas(f3, fpath)
end
