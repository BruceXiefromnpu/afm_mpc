% This script plots the LQR based root locus for the "lowgain"
% version of things.

clear, clc
saveon = true;
plants = CanonPlants.plants_ns14();
G = plants.SYS;
G_recyc = plants.sys_recyc;
Ts = G_recyc.Ts;

const_sig = false;

figbase = 30;

lqr_locus_figfile = 'lqr_locus_constsig_0p9.svg';
cmplx_rad = 0.9;
[Q1, R0, S1, P_x] = build_control_constsigma(G_recyc, cmplx_rad);

f3 = figure(3+figbase); clf

t = (0:.01:pi);
x = cmplx_rad*sin(t);
y = cmplx_rad*cos(t);
plot(real(P_x), imag(P_x), 'ok')
hold on
plot(x, y, 'k')

gam_min = 0.001;
gam_max = 10000;

ax = gca();
rgb1 = [0.230, 0.299, 0.754];
rgb2 = [0.706, 0.016, 0.150];
s_ = linspace(0,1, 400);
color_map = diverging_map(s_, rgb1, rgb2);



[ax, C_hand] = lqr_locus(G_recyc, Q1, 1, S1, gam_min, gam_max, ax,...
  'color_map', color_map);
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
  saveas(f3, fpath);
end



  








