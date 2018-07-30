clear
clc
addpath('functions')
fig_path1 = fullfile('latex/figures', 'timeopt_slowpz_illustrate_y.svg');
fig_path2 = fullfile('latex/figures', 'timeopt_slowpz_illustrate_u.svg');
load /media/labserver/mpc-journal/x-axis_sines_info_out_2-8-2018-01.mat
% load /media/labserver/mpc-journal/sysID/x-axis_sines_info_matchTs_out_2-14-2018-01.mat
whos

G_xpow = modelFit.frf.G_xpow_cors;
freqs = modelFit.frf.freq_s;



P_pow2stage = squeeze(G_xpow(1,2, :)./G_xpow(2,2, :));
P_uz2stage = squeeze(G_xpow(1,3, :)./G_xpow(3,3, :));
P_uz2pow = squeeze(G_xpow(2,3, :)./G_xpow(3,3, :));

F1 = figure(1); clf
axes('FontName','Times New Roman') % Set axis font style
box('on'); % Define box around whole figure
hold on;
width = 4.5;
height = 4;
set(F1, 'Units', 'Inches', 'Position', [0, 0, width, height],...
    'PaperUnits', 'Inches', 'PaperSize', [width, height])
set(F1, 'Color', 'w');
frfBode(P_uz2stage, freqs, F1, 'Hz', 'r');


sys = modelFit.models.G_uz2stage;
sys.InputDelay = 0;
% sys = sys*modelFit.models.G_real_extra;

frfBode(sys, freqs, F1, 'Hz', '--k')

plotPZ_freqs(sys, F1, 0, 0)
%% Create system without slow pz-pair
% grab g_drift;
Ts= sys.Ts;

pl = pole(sys);
zr = zero(sys);
pl_real = pl(find(imag(pl)==0));
zr_real = zr(find(imag(zr)==0));
gdrift = zpk(zr_real(1), pl_real(end), 1, Ts);
gdrift = gdrift/dcgain(gdrift);

sys_fast = minreal((1/gdrift)*sys)
frfBode(sys_fast, freqs, F1, 'g', 'hz')

%
%%
sys_recyc = SSTools.deltaUkSys(sys);
% sysct = tf(100^2, [1 2*0.1*100, 100^2]);
% sys = ss(c2d(sysct, 40e-6));
N = 250;
dumax = 0.3;
to_bisect = TimeOptBisect(sys_recyc, dumax);
[Nx, Nu] = SSTools.getNxNu(sys_recyc);

yref = 2;
x0 = Nx*0;
xf = Nx*yref;
[xx, uu] = to_bisect.time_opt_bisect(x0, xf, 'verbose', true);
nu = length(uu.Data);

u = [uu.Data; ones(N-nu, 1)*Nu*yref];
t = [0:1:N-1]'*Ts;
y = lsim(sys_recyc, u, t);

% without slow-pz
sys_recyc_fast = SSTools.deltaUkSys(sys_fast);
to_bisect_fast = TimeOptBisect(sys_recyc_fast, dumax);
[Nx, Nu] = SSTools.getNxNu(sys_recyc_fast);

x0 = Nx*0;
xf = Nx*yref;
[xx_fast, uu_fast] = to_bisect_fast.time_opt_bisect(x0, xf, 'verbose', true);
nuf = length(uu_fast.Data);

u_fast = [uu_fast.Data; ones(N-nuf, 1)*Nu*yref];
t = [0:1:N-1]'*Ts;
y_fast = lsim(sys_recyc_fast, u_fast, t);

%%
F4 = figure(4); 
axes('FontName','Times New Roman') % Set axis font style
box('on'); % Define box around whole figure
hold on;
width = 2.5;
height = 2.5;
set(F4, 'Units', 'Inches', 'Position', [0, 0, width, height],...
    'PaperUnits', 'Inches', 'PaperSize', [width, height])
set(F4, 'Color', 'w');



clf; hold on
plot(t, u);
plot(t, u_fast)
yl = ylabel('control $\Delta u(k)$'); %, 'interpreter', 'latex')
xl = xlabel('time [s]')

F3 = figure(3); clf
axes('FontName','Times New Roman') % Set axis font style
box('on'); % Define box around whole figure
hold on;
width = 2.5;
height = 2.5;
set(F3, 'Units', 'Inches', 'Position', [0, 0, width, height],...
    'PaperUnits', 'Inches', 'PaperSize', [width, height])
set(F3, 'Color', 'w');

plot(t, y)

plot(t, y_fast);
ylabel('output y(k)')
xlabel('time [s]')


if 1
    saveas(F3, fig_path1)
    saveas(F4, fig_path2)
end