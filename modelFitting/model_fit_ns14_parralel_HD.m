clear
clc


%%
% ------------------- Fit Hysteresis + Sat -------------------------------------

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

umax = (abs(ux));
ymax = abs(min(yx));


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



clear PIHyst
[r2, w2] = PIHyst.fit_hyst_weights(downsample(ux, 100), downsample(yprime, 100), Nhyst);
u_pi  = PIHyst.hyst_play_op(ux, r2, w2, w*0);
y_hyst = lsim(Gvib*gdrift, u_pi, tvec);
fprintf('classic hyst residual: %.2f\n', sum((y_hyst - yx).^2));


figure(1000); clf
plot(tvec, yprime)
hold on
plot(tvec, u_pisat)
grid on
%%
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
% theta0 = [w; ws]%+.05
% clc
% nsd = 4; % drift model size.
% gd = canon(gdrift, 'modal');
% gd_ = gd;
% Gvib_ =  ss(Gvib);
% % gd_ = d2d(gd, gd.Ts*10);
% % Gvib_ = ss(d2d(Gvib, Gvib.Ts*1));
% ux_ = downsample(ux, 100);
% yprime_ = downsample(yx, 1);
% % gd.d 
% gd_.c = gd_.c.*gd_.b';
% gd_.b = gd_.b*0 + 1;
% 
% lams = [diag(gd_.a); 0.8; 0.9];
% theta_gd0 = [lams; [gd_.c'; 0.0001; 0.0001]; gd_.d];
% theta_hyst0 = [w];
% theta_sat0 = [ws];
% theta0 = [theta_hyst0; theta_sat0; theta_gd0];
% lb = repmat([-Inf], 1, Nhyst+Nsat+2*nsd+1);
% lb(end-2*nsd:end-2*nsd+nsd-1) = -1+1e-9;
% ub = -lb;
% [lb', theta0]
% lb(1:Nhyst) = 0;
% lb(end) = 0;
nsd = 4; % drift model size.
gd = canon(gdrift, 'modal');
% gd.d 
gd.c = gd.c.*gd.b';
gd.b = gd.b*0 + 1;

theta_gd0 = [[diag(gd.a); 0.998; 0.99]; [gd.c'; 0.0001; 0.0001]; gd.d*0.5];
theta_hyst0 = [w];
theta_sat0 = [ws];
theta0 = [theta_hyst0; theta_sat0; theta_gd0];

lb = repmat([-Inf], 1, Nhyst+Nsat+2*nsd+1);
lb(end-2*nsd:end-2*nsd+nsd-1) = -1+1e-9;
ub = -lb;
[lb', theta0]


Jfunc = @(theta) hyst_drift_paralel_der(theta, ux, yx, tvec, Nhyst, Nsat, nsd, absorbDelay(Gvib), r, d);

HSopts = optimoptions(@lsqnonlin, 'Display', 'iter',...
    'FunctionTolerance', 1e-7, 'MaxIter', 5000,'MaxFunctionEvaluations', 5000,...
    'StepTolerance', 1e-8, 'Jacobian','on', 'CheckGradients', false);
  
theta = lsqnonlin(Jfunc, theta0, lb, ub, HSopts);
%%
w_hyst = theta(1:Nhyst);
w_sat = theta(Nhyst+1:Nhyst+Nsat);
theta_gd = theta(Nhyst+1+Nsat:end);
a = diag(theta_gd(1:nsd));
c = theta_gd(nsd+1:end-1)';
b = c'*0+1;
gd_fit = ss(a, b, c, theta_gd(end), Gvib.Ts);

u_drft = lsim(gd_fit, ux, tvec);

u_hyst = PIHyst.hyst_play_sat_op(ux, r, w_hyst, d, w_sat, w_hyst*0);
u_effec = u_drft+u_hyst;
y_fit_mine = lsim(Gvib, u_effec, tvec);

figure(2000); 
plot(tvec, y_fit_mine, '--k')
hold on
plot(tvec, yx)
legend('y-exp', 'y-hyst-sat', 'y-hyst', 'parallel')

% figure(1000); clf
% plot(tvec, u_effec, 'm')
% hold on


%%
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
%%
% c
[theta, theta0]
[~, y_fit] =  err_func(theta);

figure(2000)
plot(tvec, y_fit, '--b')
legend('y-exp', 'y-hyst-sat', 'y-hyst', 'hyst+sat + drift parallel')

figure(3000)
plot(tvec, yx - y_fit, '--b')
legend('y-exp', 'y-hyst-sat', 'y-hyst', 'hyst+sat + drift parallel')

fprintf('-------------------------------------------------------------\n')
fprintf('results summary\n')
fprintf('classic hyst+sat residual: %.2f\n', sum((y_hyst_sat - yx).^2));
fprintf('classic hyst residual: %.2f\n', sum((y_hyst - yx).^2));
fprintf('new parallel hyst+sat residual: %.2f\n', sum((y_fit - yx).^2));



w_hyst = theta(1:Nhyst);
w_sat = theta(Nhyst+1:Nhyst+Nsat);
theta_gd = theta(Nhyst+1+Nsat:end);

a = diag(theta_gd(1:nsd));
c = theta_gd(nsd+1:end-1)';
b = c'*0+1;

gd = ss(a, b, c, theta_gd(end), Gvib.Ts);

hyst_drift_paral = struct('gdrift', gd, 'w', w_hyst, 'r', r,...
                         'ws', w_sat, 'd', d);

modelFit.models.hyst_drift_paral = hyst_drift_paral;
% % % % modelFit.models.G_uz2powI = G_deluz2Ipow*g_der;
% % % % modelFit.models.G_deluz2powI = G_deluz2Ipow;
% % % % modelFit.models.g_deluz2pow_1norm = nm1;
% % % % modelFit.models.du_max_nm1 = delumax;
% % % modelFit.models.G_uz2stage = sys_stage_log;
% % % modelFit.models.G_uz2powI = G_deluz2Ipow*g_der;
% % % modelFit.models.G_deluz2powI = G_deluz2Ipow;
% % % modelFit.models.g_deluz2pow_1norm = nm1;
% % % modelFit.models.du_max_nm1 = delumax;
% % % modelFit.models.Gvib = Gvib;
% % % modelFit.models.gdrift_1p0 = gdrift;
% % % modelFit.models.hyst = hyst;
% % % modelFit.models.hyst_sat = hyst_sat;
if 1
    save(modelFit_file, 'modelFit');
end




