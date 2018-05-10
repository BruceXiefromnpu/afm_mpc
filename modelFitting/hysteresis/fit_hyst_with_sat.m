% This script tries to fit ONly the hysteresis model. It uses the separeted
% Gvib and drift model fit in fit_drift.m

clear

load('hystID_data_5-4-2018_01.mat')

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



%%
% start with prior infomation for guess
guess_data = load('steps_hyst_model.mat');
if 1
  r = guess_data.r;
  w = guess_data.theta(1:7);
  d = guess_data.d;
  ws = guess_data.theta(8:end);
end
umax = max(abs(ux));
nw = 7;

% r = linspace(0, umax*0.5, nw);
% r = ([0:nw-1]'./(nw) )*umax
% d = [0; 4.5; 9];

% d = [0; 5; 10]
% w = ones(nw, 1);

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


%%
figure(18); clf;
h1 = plot(downsample(tvec, 50), downsample(yx, 50));
hold on
h2 = plot(downsample(tvec, 50), downsample(y_est_all, 50), 'r');
h1.DisplayName = 'Experimental';
h2.DisplayName = 'Fit';
legend([h1, h2])
grid on

w = theta(1:nw);
ws = theta(nw+1:end);
[rp, wp, dp, wsp] = PIHyst.invert_hyst_sat_PI(r, w, d, ws);



%%
save('steps_hyst_model_withsat.mat', 'rp', 'wp', 'r', 'd', 'dp', 'wsp',...
  'w', 'ws', 'theta', 'umax', 'Gvib', 'gdrift')



%%
% ----------------------------------------------------------------------- %
% Now, we try fitting the drift and hysteresis together, without this
% saturation business
% w = ones(nw, 1);
dat = load('steps_hyst_model.mat')
gdrift = modelFit.models.drift.gdrift;
w0 = [3.9355, 4.9039, -4.701, 2.7965, 1.2377, -2.8405, 2.5493];
idx_w0 = [1:length(w0)];
ws0 = ws(:)';
idx_ws0 = [1:length(ws0)] + idx_w0(1);
[num, den] = tfdata(gdrift, 'v');

idx_num = [1,2,3] + idx_ws0(end);
idx_den = [1,2] + idx_num(end);

theta_0 = [w0, ws0, num, den(2:end)];

%%
clc
gdrift_h = @(theta) tf(theta(idx_num), [1, theta(idx_den)], Ts);
G = @(theta) Gvib*gdrift_h(theta);
u_hyst_fun = @(theta) PIHyst.hyst_play_sat_op(ux, r, theta(idx_w0), d, theta(idx_ws0), w0*0);

cost_hyst_only = @(theta) downsample(lsim(G(theta), u_hyst_fun(theta),... 
                tvec) - yx, 50); 
              
              
opts = optimoptions(@lsqnonlin);
opts.MaxFunctionEvaluations = 10000;
opts.Display = 'iter';

lb = [];
ub = [];
theta = lsqnonlin(cost_hyst_only, theta_0, lb, ub, opts);

% theta = theta_0;
%%
u_est_steps_all = PIHyst.hyst_play_sat_op(ux, r, theta(idx_w0), d, theta(idx_ws0), w*0);

y_est_all = lsim(G(theta), u_est_steps_all, tvec);



figure(19); clf;
h1 = plot(downsample(tvec, 50), downsample(yx, 50));
hold on
h2 = plot(downsample(tvec, 50), downsample(y_est_all, 50), 'r');
h1.DisplayName = 'Experimental';
h2.DisplayName = 'Fit';
legend([h1, h2])
grid on

w = theta(1:nw);
ws = theta(nw+1:end);
[rp, wp, dp, wsp] = PIHyst.invert_hyst_sat_PI(r, w, d, ws);



%%
save('steps_hyst_model.mat', 'rp', 'wp', 'r', 'd', 'dp', 'wsp',...
  'w', 'ws', 'theta', 'umax', 'Gvib', 'gdrift')


