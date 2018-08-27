

addpath(fullfile(getMatPath(), 'afm_mpc_journal', 'functions'))
addpath(fullfile(getMatPath(), 'afm_mpc_journal', 'functions', 'canon'))

clear, clc
saveon = true;
plants = CanonPlants.plants_ns14();
G = plants.SYS;
G2 = plants.Gvib;


p1 = pole(G);
z1 = tzero(G2);
p2 = pole(absorbDelay(G2));

assert( all(abs(sort(p1) - sort(p2)) < 1e-13) == true)
%%
fig = mkfig(1, 1.75, 2); clf
%                                           lft   bt   wd   ht
ax = axes('Units', 'inches', 'Position', [0.33,  0.2, 1.4 1.7500]);
ax = gca();
hold on
plot(real(p2), imag(p2), 'xb')
plot(real(z1), imag(z1), 'ob')

title('')
xlim([0.75, 1.001])
ylim([0, 0.3])

zgrid([.1, .3, .5, .7], [.1, .25, .5])
hold on

xlab = xlabel('Re', 'FontSize', 8, 'Interpreter', 'latex');
ylab = ylabel('Im', 'FontSize', 8, 'Interpreter', 'latex');
set(ylab, 'Units', 'inches', 'Position', [-.2, 1.5/2, 0]);
set(xlab, 'Units', 'inches', 'Position', [1.5/2, -.1, 0]);

saveas(ax, fullfile(PATHS.jfig, 'pzplot.svg'))
% tighten_axis(fig, ax)
%%
S = ss2tex(G)

fname = fullfile(PATHS.MPCJ_root, 'latex', 'Gvib_data.tex');
fid = fopen(fname, 'w+');
fprintf(fid, '%s', S);

fclose(fid);