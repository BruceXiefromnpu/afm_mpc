% This script tries to fit ONly the hysteresis model. It uses the separeted
% Gvib and drift model fit in fit_drift.m

clear

load('hystID_data_4-30-2018_01.mat')

% To remove the vibrational and drift aspects, we need those models. Load them:
modelFit_file = fullfile(PATHS.sysid, 'FRF_data_current_stage2.mat');
load(modelFit_file)
Ts = modelFit.frf.Ts;
Gvib = modelFit.models.drift.Gvib;
gdrift = modelFit.models.drift.gdrift;

ms = 1e3;

ux = hystData.u_exp;
yx = hystData.y_exp - mean(hystData.y_exp(1:500));
tvec = hystData.t_exp;

% ux = downsample(ux, 50);
figure
plot(ux)

% %%
% 
% umax = max(abs(ux));
% n = 7;
% % r = linspace(0, umax, n);
% r = ([0:n-1]'./(n) )*umax/2
% w = ones(n, 1);
% 
% % hyst_fun = @(w, y0) PIHyst.hyst_play_op(u_exp(idx_ureset_end+1:end), r, w, y0);
% 
% G = Gvib*gdrift;
% 
% % u_hyst_fun = @(theta) PIHyst.hyst_play_op(ux, theta(n+1:end), theta(1:n), w*0);
% u_hyst_fun = @(theta) PIHyst.hyst_play_op(ux, r, theta(1:n), w*0);
% 
% cost_hyst_only = @(theta) downsample(lsim(G, u_hyst_fun(theta),... 
%                 tvec) - yx, 50); 
%               
%               
% opts = optimoptions(@lsqnonlin);
% opts.MaxFunctionEvaluations = 10000;
% opts.Display = 'iter';
% % theta_0 = [w; r];
% theta_0 = w;
% lb = [0*w];
% ub = [];
% 
% theta = lsqnonlin(cost_hyst_only, theta_0, lb, ub, opts);
% 
% % u_est_steps_all = PIHyst.hyst_play_op(ux, theta(n+1:end), theta(1:n), w*0);
% u_est_steps_all = PIHyst.hyst_play_op(ux, r, theta(1:n), w*0);
% 
% y_est_all = lsim(G, u_est_steps_all, tvec);
%%
umax = max(abs(ux));
nw = 7;
% r = linspace(0, umax, n);
r = ([0:nw-1]'./(nw) )*umax
d = [0; 4.5; 9];
w = ones(nw, 1);
% ws = 0.5*ones(2,1);
ws = [1; 0.0001; 0.0001];
G = Gvib*gdrift;
u_hyst_fun = @(theta) PIHyst.hyst_play_sat_op(ux, r, theta(1:nw), d, theta(nw+1:end), w*0);

cost_hyst_only = @(theta) downsample(lsim(G, u_hyst_fun(theta),... 
                tvec) - yx, 50); 
              
              
opts = optimoptions(@lsqnonlin);
opts.MaxFunctionEvaluations = 10000;
opts.Display = 'iter';
theta_0 = [w; ws];
% theta_0 = w;
% lb = [0*w; 0*ws];
lb = [];
ub = [];
theta = lsqnonlin(cost_hyst_only, theta_0, lb, ub, opts);

% theta = theta_0;
u_est_steps_all = PIHyst.hyst_play_sat_op(ux, r, theta(1:nw), d, theta(nw+1:end), w*0);

y_est_all = lsim(G, u_est_steps_all, tvec);


%
%%
figure(18); clf;
h1 = plot(downsample(tvec, 50), downsample(yx, 50));
hold on
h2 = plot(downsample(tvec, 50), downsample(y_est_all, 50), 'r');
% h3 = plot(tvec, ux*dcgain(G), 'LineWidth', 1.5);

h1.DisplayName = 'Experimental';
h2.DisplayName = 'Fit';
% h3.DisplayName = 'u*dcgain(G)';
% legend([h1, h2, h3])
legend([h1, h2])
grid on




%%
[rp, wp] = PIHyst.invert_hyst_PI(r, theta);

save('steps_hyst_model.mat', 'rp', 'wp', 'r', 'theta_hyst', 'umax', 'Gvib', 'gdrift')




