clear
clc
ms = 1e3;
modelFit_file = fullfile(PATHS.sysid, 'FRF_data_current_stage2.mat');
load(modelFit_file)
load('driftID_data_4-30-2018_01.mat')
Ts = modelFit.frf.Ts;
G_uz2stage = modelFit.models.G_uz2stage;



k_start = 2677;
y_exp = driftData.y_exp(k_start:end) - mean(driftData.y_exp(1:k_start));
u_exp = driftData.u_exp(k_start:end);
t_vec = (0:length(y_exp)-1)'*Ts;

figure(100); clf
plot_ind(t_vec, y_exp)

theta0 = [0.9922    0.9997    0.9997    0.9927    .8];
lb = [-1, -1, -1, -1, -Inf];
ub = -lb;
np = 2;

Gvib = eject_gdrift(G_uz2stage);
gdrift_cost = @(theta)fit_gdrift(theta, Gvib, y_exp, u_exp, t_vec, np);

opts = optimoptions(@lsqnonlin);
opts.MaxFunctionEvaluations = 5000;
opts.MaxIterations = 1000;
opts.Display = 'iter';
theta = lsqnonlin(gdrift_cost, theta0, lb, ub, opts);

gdrift = zpk(theta(np+1:end-1), theta(1:np), theta(end), Ts);

ydrift_est0 = lsim(Gvib*gdrift, u_exp, t_vec);
ynom = lsim(G_uz2stage, u_exp, t_vec);
clf
plot(t_vec, y_exp)
hold on, grid on
plot(t_vec, ydrift_est0, '--r')


%%
modelFit.models.drift.Gvib = Gvib;
modelFit.models.drift.gdrift = gdrift;

save(modelFit_file, 'modelFit')







