clear
clc
addpath('functions')

fig_path = fullfile('latex/figures', 'frf_xdir.svg');
load /media/labserver/mpc-journal/x-axis_sines_info_out_2-8-2018-01.mat
% load /media/labserver/mpc-journal/sysID/x-axis_sines_info_matchTs_out_2-14-2018-01.mat
whos
saveon = 0;
G_xpow = modelFit.frf.G_xpow_cors;
freqs = modelFit.frf.freq_s;



P_pow2stage = squeeze(G_xpow(1,2, :)./G_xpow(2,2, :));
P_uz2stage = squeeze(G_xpow(1,3, :)./G_xpow(3,3, :));
P_uz2pow = squeeze(G_xpow(2,3, :)./G_xpow(3,3, :));

F1 = figure(1); clf
axes('FontName','Times New Roman'); % Set axis font style
box('on'); % Define box around whole figure
hold on;
width = 4.5;
height = 4;
set(F1, 'Units', 'Inches', 'Position', [0, 0, width, height],...
    'PaperUnits', 'Inches', 'PaperSize', [width, height])
set(F1, 'Color', 'w');
frfBode(P_uz2stage, freqs, F1, 'r', 'Hz');


sys = modelFit.models.G_uz2stage;
sys.InputDelay = 9;
sys = sys*modelFit.models.G_real_extra;

frfBode(sys, freqs, F1, '--k', 'Hz');

plotPZ_freqs(sys, F1, 0, 0);

if saveon
    saveas(F1, fig_path)
end
%%
clc
F2 = mkfig(2, 3.5, 3);
G_uz2pow = modelFit.models.G_uz2pow;
frfBode(P_pow2stage, freqs, F2, 'r', 'Hz');
frfBode(G_uz2pow, freqs, F2, 'b', 'Hz');

subplot(2,1,1)
ylim([-80, 30])
ax1 = gca;
set(ax1, 'YMinorTick', 'on', 'XMinorGrid', 'on', 'XMinorTick', 'on')
set(ax1, 'XTickLabel', {'$10^1$', '$10^2$', '$10^3', '$10^4$'}, 'XTick',...
    [10, 100, 1000, 10000])

subplot(2,1,2)
ax2 = gca;
ylim([-1000, 50])
set(ax2, 'XTickLabel', {'$10^1$', '$10^2$', '$10^3', '$10^4$'}, 'XTick',...
    [10, 100, 1000, 10000])
set(ax2, 'YMinorTick', 'on', 'XMinorGrid', 'on', 'XMinorTick', 'on')
if saveon
    saveas(F2, 'latex/figures/G_uz2pow.svg')
end