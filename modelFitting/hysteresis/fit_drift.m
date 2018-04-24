clear
clc

modelFit_file = fullfile(PATHS.sysid, 'FRF_data_current_stage2.mat');
load(modelFit_file)

Ts = modelFit.frf.Ts;
G_uz2stage = modelFit.models.G_uz2stage;


ms = 1e3;

%Case 1:

clc
To_half = 28*Ts;
To = 2*To_half;
wo = 2*pi/To;

% Build the u-reset.
t1 = 18;
tf = 100;
umax = 5;
k1 = 0.15;
u_reset = PIHyst.gen_reset_u(t1, tf, Ts, k1, umax);
figure(100)
plot(u_reset)
%%

slewfname_in = 'drift_id_datain_4-23-2018_02.csv';
slewfpath_in = fullfile(PATHS.exp, 'sysID', slewfname_in);

slewfname_out = 'drift_id_4-23-2018_02.csv';
slewfpath_out = fullfile(PATHS.exp, 'sysID', slewfname_out);

dat = csvread(slewfpath_out);
%%
idx_reset = length(u_reset);
k_start = idx_reset+474068 + 19959;
u_exp = dat(k_start:end, 1);
y_exp = dat(k_start:end, 2);


k2 = find(diff(u_exp) ~=0, 1, 'first');
y_exp = y_exp - mean(y_exp(1:k2));
k3 = 13567;
y_exp = y_exp(1:k3);
u_exp = u_exp(1:k3);
t_vec = (0:length(y_exp)-1)'*Ts;

figure(100); clf
plot(t_vec, y_exp)
%%
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
%%
gdrift = zpk(theta(np+1:end-1), theta(1:np), theta(end), Ts);

ydrift_est0 = lsim(Gvib*gdrift, u_exp, t_vec);
ynom = lsim(G_uz2stage, u_exp, t_vec);
clf
plot(t_vec, y_exp)
hold on, grid on
plot(t_vec, ydrift_est0, '--r')


%%

save('drift_data.mat', 'G_uz2stage', 'Gvib', 'gdrift')







