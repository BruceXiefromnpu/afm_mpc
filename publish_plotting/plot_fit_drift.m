% --------------------- Now, Fit the drift model -----------------
clear
addpath(fullfile(getMatPath, 'afm_mpc_journal', 'functions'))
addpath(fullfile(getMatPath, 'afm_mpc_journal', 'functions', 'canon'))

% load(fullfile(PATHS.sysid, 'hysteresis', 'driftID_data_06-05-2018_01_amp_1p0.mat'))
% load(fullfile(PATHS.sysid, 'hysteresis', 'driftID_data_5-29-2018_01.mat'))
saveon = true;
[plants, frf_data, MF] = CanonPlants.plants_ns14(9, 2);
fpath_driftID = MF.heritage.fpath_driftID;

% make sure the data hash not changed
drift_hash_expected = MF.heritage.drift_data_hash;
drift_hash_current = DataHash(fpath_driftID, MF.heritage.hash_opts);
assert(strcmp(drift_hash_expected, drift_hash_current))

load(fpath_driftID)

% fullfile(PATHS.sysid, 'hysteresis', 'driftID_data_9-10-2018_01.mat'))

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

h0 = plot(ax, t_exp, u_exp, '-.k');
h0.DisplayName = 'Input';

h1 = plot(t_exp, y_exp, '-b');
h1.DisplayName = 'Exp. Step Response';

h2 = plot(t_exp, ydrift_est0, '--r');
h2.DisplayName = '$\hat{G}_{\textrm{vib}}\hat{G}_d$';

h3 = plot(t_exp, y_vib, ':k');
h3.DisplayName = '$\hat{G}_{\textrm{vib}}$';
leg1 = legend([h0, h1, h2, h3]);
xlim([0, 0.28])
ylim([-0.005, 1.05])
% set(ax, 'YTick', (0:0.025:0.15))

grid on
% grid minor;
xlabel('time [s]')
ylabel('$y_X$ [v]')
ax = gca;

tighten_axis(fig, ax)
set(leg1, 'Units', 'inches', 'Box', 'off', 'Position', [1.7751 0.5873 1.5782 0.6273]);
if saveon
  saveas(fig, fullfile(PATHS.jfig(), 'drift_fit.svg'))
end