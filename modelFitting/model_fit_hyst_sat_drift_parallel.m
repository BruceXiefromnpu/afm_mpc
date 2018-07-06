clear
clc


%%
% ------------------- Fit Hysteresis + Sat -------------------------------------
addpath(fullfile(getMatPath(), 'afm_mpc_journal', 'functions'));
modelFit_file = fullfile(PATHS.sysid, 'FRF_data', 'x-axis_sines_infoFourierCoef_5-30-2018-01.mat');

load(modelFit_file)
Gvib = modelFit.models.Gvib;
gdrift = modelFit.models.gdrift;

fprintf('============================================\n')
hyst_file = 'hystID_data_5-4-2018_01.mat';
% hyst_file = 'hystID_data_6-1-2018_01.mat';

% hyst_file = 'hystID_data_6-1-2018_01.mat';
hyst_path = fullfile(PATHS.sysid, 'hysteresis', hyst_file);
load(hyst_path)

kk = length(hystData.y_exp);
ux = hystData.u_exp(1:kk);
yx = hystData.y_exp(1:kk) - mean(hystData.y_exp(1:500));
tvec = hystData.t_exp(1:kk);

umax = max(abs(ux));
ymax = max(abs(yx));


Nhyst = 7;
nw = Nhyst;
Nsat = 7;

yprime = lsim(1/(gdrift*dcgain(Gvib)), yx, tvec);
[r, w, d, ws] = PIHyst.fit_hyst_sat_weights(downsample(ux, 100), downsample(yprime, 100), Nhyst, Nsat);
u_pisat = PIHyst.hyst_play_sat_op(ux, r, w, d, ws, w*0);
y_hyst_sat = lsim(Gvib*gdrift, u_pisat, tvec);
fprintf('classic hyst+sat residual: %.2f\n', sum((y_hyst_sat - yx).^2));

[rp, wp] = PIHyst.invert_hyst_PI(r, w);
[dp, wsp] = PIHyst.invert_sat(d, ws);


hyst_sat = struct('r', r, 'w', w, 'rp', rp, 'wp', wp,...
  'd', d, 'ws', ws, 'dp', dp, 'wsp', wsp);


[r2, w2] = PIHyst.fit_hyst_weights(downsample(ux, 100), downsample(yprime, 100), Nhyst);
u_pi  = PIHyst.hyst_play_op(ux, r2, w2, w*0);
y_hyst = lsim(Gvib*gdrift, u_pi, tvec);
fprintf('classic hyst residual: %.2f\n', sum((y_hyst - yx).^2));


figure(1000); clf
plot(tvec, yprime)
hold on
plot(tvec, u_pisat)
grid on

figure(2000); clf
plot(tvec, yx)
hold on
plot(tvec, y_hyst_sat)
plot(tvec, y_hyst, '--')

legend('y-exp', 'y-hyst-sat', 'y-hyst')
grid on

figure(3000); clf
plot(tvec, y_hyst_sat - yx)
hold on
plot(tvec, y_hyst-yx, '--')
legend('y-hyst-sat error', 'y-hyst error')
title('model error')
grid on

[rp2, wp2] = PIHyst.invert_hyst_PI(r2, w2);
hyst = struct('r', r2, 'w', w2, 'rp', rp2, 'wp', wp2);
%%
% ---------- Fit Hyst[] * Sat [] + drift (parallel structure) -------------

% construct the intervals
n_d = (Nsat - 1)/2;
r = ([0:Nhyst-1]'./(Nhyst) )*umax;
% id_plus = (1:n_d);
% id_neg = (-n_d:-1);
% dplus = ((id_plus - 0.5)/n_d ) * ymax;
% dmin = ( (id_neg + 0.5)/n_d ) *ymax; 
% d = [dmin, 0, dplus]';


gd = canon(gdrift, 'modal');
gd.c = gd.c.*gd.b'; % scale
gd.b = gd.b*0 + 1;

has_feedthrough = false;

if has_feedthrough
  a_ = [diag(gd.a); 0.998; 0.99];
  nsd = length(a_);
  theta_gd0 = [a_; [gd.c'; 0.0001; 0.0001]; gd.d*0.5];
  theta_hyst0 = [w];
  theta_sat0 = [ws];
  theta0 = [theta_hyst0; theta_sat0; theta_gd0];
  
  % Lower and upper bounds to enforce stability of gdrift.
  lb = repmat([-Inf], 1, Nhyst+Nsat+2*nsd+1);
  lb(end-2*nsd:end-2*nsd+nsd-1) = -1+1e-9;
  ub = -lb;
  [lb', theta0]
else
  % a_ = [diag(gd.a); 0.998; 0.99];
  % c_ = [gd.c'; 0.0001; 0.0001]
  a_ = [diag(gd.a)];
  c_ = gd.c';
  nsd = length(a_);
  theta_gd0 = [a_; c_];
  theta_hyst0 = [w]*gd.d;
  theta_sat0 = [ws];
  theta0 = [theta_hyst0; theta_sat0; theta_gd0];
  
  % Lower and upper bounds to enforce stability of gdrift.
  lb = repmat([-Inf], 1, Nhyst+Nsat+2*nsd);
  lb(Nsat+Nhyst+1:Nsat+Nhyst+nsd) = -1+1e-9;
  ub = -lb;
  [lb', theta0]  
  
  
end
%
Jfunc = @(theta) hyst_drift_paralel_der(theta, ux, yx, tvec, Nhyst, Nsat, nsd, absorbDelay(Gvib), r, d, has_feedthrough);

HSopts = optimoptions(@lsqnonlin, 'Display', 'iter',...
    'MaxIter', 5000,'MaxFunctionEvaluations', 5000,...
    'Jacobian','on', 'CheckGradients', false);
  
theta = lsqnonlin(Jfunc, theta0, lb, ub, HSopts);
%%
w_hyst = theta(1:Nhyst);
w_sat = theta(Nhyst+1:Nhyst+Nsat);
theta_gd = theta(Nhyst+1+Nsat:end);
a = diag(theta_gd(1:nsd));
if has_feedthrough
  c = theta_gd(nsd+1:end-1)';
  b = c'*0+1;
  D = theta_gd(end);
else
  c = theta_gd(nsd+1:end)';
  b = c'*0+1;
  D = 0;
end 
%%
gd_fit = ss(a, b, c, D, Gvib.Ts);

[rp, wp] = PIHyst.invert_hyst_PI(r, w_hyst);
[dp, wsp] = PIHyst.invert_sat(d, w_sat);


u_drft = lsim(gd_fit, ux, tvec);

u_hyst = PIHyst.hyst_play_sat_op(ux, r, w_hyst, d, w_sat, w_hyst*0);
u_effec = u_drft+u_hyst;
y_fit_mine = lsim(Gvib, u_effec, tvec);


figure(2000)
plot(tvec, y_fit_mine, '--b')
legend('y-exp', 'y-hyst-sat', 'y-hyst', 'hyst+sat + drift parallel')

figure(3000)
plot(tvec, yx - y_fit_mine, '--m')
legend('y-exp', 'y-hyst-sat', 'y-hyst', 'hyst+sat + drift parallel')

fprintf('----------------------------------------------\n')
fprintf('------------results summary-------------------\n')
fprintf('classic hyst+sat residual: %.2f\n', sum((y_hyst_sat - yx).^2));
fprintf('classic hyst residual: %.2f\n', sum((y_hyst - yx).^2));
fprintf('new parallel hyst+sat residual: %.2f\n', sum((y_fit_mine - yx).^2));


hyst_drift_paral = struct('gdrift', gd_fit, 'w', w_hyst, 'r', r, 'rp', rp,...
                         'wp', wp, 'ws', w_sat, 'd', d, 'dp', dp, 'wsp', wsp);

modelFit.models.hyst_drift_paral = hyst_drift_paral;

if 1
    save(modelFit_file, 'modelFit');
end

%% ------------------------------------------------------------------------
% Old method to fit the paralel model without the derivatives.

% nsd = 4; % drift model size.
% gd = canon(gdrift, 'modal');
% % gd.d 
% gd.c = gd.c.*gd.b';
% gd.b = gd.b*0 + 1;
% 
% theta_gd0 = [[diag(gd.a); 0.998; 0.99]; [gd.c'; 0.0001; 0.0001]; gd.d*0.5];
% theta_hyst0 = [w];
% theta_sat0 = [ws];
% theta0 = [theta_hyst0; theta_sat0; theta_gd0];
% 
% 
% lb = repmat([-Inf], 1, Nhyst+Nsat+2*nsd+1);
% 
% lb(end-2*nsd:end-2*nsd+nsd-1) = -1;
% [lb', theta0]
clc
err_func = @(theta) hyst_drift_paralel_der(theta, ux, yx, tvec, Nhyst, Nsat, nsd, absorbDelay(Gvib), r, d);

err_func(theta0);
opts = optimoptions(@lsqnonlin);
opts.Display = 'iter';

theta = lsqnonlin(err_func, theta0, lb, -lb, opts)

[theta, theta0]
[~, y_fit] =  err_func(theta);



