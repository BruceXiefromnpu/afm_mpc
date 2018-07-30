% --------------------- Now, Fit the drift model -----------------
clear
% load(fullfile(PATHS.sysid, 'hysteresis', 'driftID_data_06-05-2018_01_amp_1p0.mat'))
load(fullfile(PATHS.sysid, 'hysteresis', 'driftID_data_5-29-2018_01.mat'))
[plants, frf_data] = CanonPlants.plants_ns14();

Gvib = plants.SYS; %modelFit.models.G_uz2stage;
Ts = Gvib.Ts;

% G_drift = plants.gdrift;

k_start = 2677;
y_exp = driftData.y_exp(k_start:end) - mean(driftData.y_exp(1:k_start));
u_exp = driftData.u_exp(k_start:end);
t_exp = (0:length(y_exp)-1)'*Ts;

theta0 = [0.9922    0.9997    0.9997    0.9927    .8];
lb = [-1, -1, -1, -1, -Inf];
ub = -lb;
np = 2;
% normalize_dc = true;
% Gvib = eject_gdrift(G_uz2stage, normalize_dc);

gdrift_cost = @(theta)fit_gdrift(theta, Gvib, y_exp, u_exp, t_exp, np);

opts = optimoptions(@lsqnonlin);
opts.MaxFunctionEvaluations = 5000;
opts.MaxIterations = 1000;
opts.Display = 'iter';
theta = lsqnonlin(gdrift_cost, theta0, lb, ub, opts);

gdrift = zpk(theta(np+1:end-1), theta(1:np), theta(end), Ts);


ydrift_est0 = lsim(Gvib*gdrift, u_exp, t_exp);
y_vib = lsim(Gvib, u_exp, t_exp);

width = 3.5;
height = 2.5;
fig = mkfig(105, width, height); clf
ax = gca();
hold on

h0 = plot(ax, t_exp, u_exp, '-.b');
h0.DisplayName = 'Input';

h1 = plot(t_exp, y_exp, '-k');
h1.DisplayName = 'Exp. Step Response';

h2 = plot(t_exp, ydrift_est0, '--r');
h2.DisplayName = '$G_{vib}G_d$';

h3 = plot(t_exp, y_vib, ':k');
h3.DisplayName = '$G_{vib}$';
leg1 = legend([h0, h1, h2, h3]);
xlim([0, 0.28])
ylim([-0.005, 0.17])
% set(ax, 'XTick', (0:0.05:0.3), 'YTick', (0:0.025:0.15))
set(ax, 'YTick', (0:0.025:0.15))

grid on
% grid minor;
xlabel('time [s]')
ylabel('$y_X$ [v]')
ax = gca;

tighten_axis(fig, ax)
set(leg1, 'Units', 'inches', 'Box', 'off', 'Position', [1.7751 0.5873 1.5782 0.6273]);

saveas(fig, fullfile(PATHS.jfig(), 'drift_fit.svg'))