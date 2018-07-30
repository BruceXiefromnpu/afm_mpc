clear
clc
ms = 1e3;
modelFit_file = fullfile(PATHS.sysid, 'FRF_data_current_stage2.mat');
load(modelFit_file)
load(fullfile(PATHS.sysid, 'hysteresis', 'driftID_data_4-30-2018_01.mat'))
Ts = modelFit.frf.Ts;
G_uz2stage = modelFit.models.G_uz2stage;



k_start = 2677;
y_exp = driftData.y_exp(k_start:end) - mean(driftData.y_exp(1:k_start));
u_exp = driftData.u_exp(k_start:end);
t_vec = (0:length(y_exp)-1)'*Ts;

F100 = figure(100); clf
plot_ind(t_vec, y_exp);


theta0 = [0.9922    0.9997    0.9997    0.9927    .8];
lb = [-1, -1, -1, -1, -Inf];
ub = -lb;
np = 2;
normalize_dc = false;
Gvib = eject_gdrift(G_uz2stage, normalize_dc);

gdrift_cost = @(theta)fit_gdrift(theta, Gvib, y_exp, u_exp, t_vec, np);

opts = optimoptions(@lsqnonlin);
opts.MaxFunctionEvaluations = 5000;
opts.MaxIterations = 1000;
opts.Display = 'iter';
theta = lsqnonlin(gdrift_cost, theta0, lb, ub, opts);

gdrift = zpk(theta(np+1:end-1), theta(1:np), theta(end), Ts);

ydrift_est0 = lsim(Gvib*gdrift, u_exp, t_vec);
y_vib = lsim(Gvib, u_exp, t_vec);

clf;
h1 = plot(t_vec, y_exp);
h1.DisplayName = 'Exp. Step Response';

hold on
h2 = plot(t_vec, ydrift_est0, '--r');
h2.DisplayName = '$G_{vib}G_d$';

h3 = plot(t_vec, y_vib, ':k');
h3.DisplayName = '$G_{vib}$';
leg1 = legend([h1, h2, h3]);
xlim([0, 0.28])
ylim([-0.005, 0.15])
grid on
% grid minor;
xlabel('time [s]')
ylabel('$y_X$ [v]')
ax = gca;
set(ax, 'XTick', (0:0.05:0.3), 'YTick', (0:0.025:0.15))
saveas(F100, fullfile(PATHS.jfig, 'drift_fit.svg'));

%%
modelFit.models.drift.Gvib = Gvib;
modelFit.models.drift.gdrift = gdrift;

save(modelFit_file, 'modelFit')







