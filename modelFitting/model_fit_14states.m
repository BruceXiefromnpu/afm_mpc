clear
clc
% The goal of this script is to get a better fit at the first complex
% resonance. It seems we need 14 states to do that.

addpath('~/gradschool/sysID/matlab/functions')

% modelFit_file = fullfile(PATHS.sysid, 'FRF_data_current_stage2.mat');
modelFit_file = fullfile(PATHS.sysid,'FRF_data', 'x-axis_sines_info_HIRESFourierCoef_5-24-2018-01.mat');

load(modelFit_file)

freqs = modelFit.frf.freq_s(:);

% G_uz2powI_frf = modelFit.frf.G_uz2powI;
G_uz2stage_frf = modelFit.frf.G_uz2stage;
G_uz2powI_frf = modelFit.frf.G_uz2powI;

omegas = freqs*2*pi;
Ts = modelFit.frf.Ts;


% --------------------------------------------------------------- %
% --------- First, we work on the stage system ----------------- %
% 
F3 = figure(3); clf
F4 = figure(4); clf
frfBode(G_uz2stage_frf, freqs, F3, 'r', 'Hz');
frfBode(G_uz2stage_frf, freqs, F4, 'r', 'Hz');

Nd2 =10;
ns2 = 14;
k_estmax = 250;
ss_opts = frf2ss_opts('Ts', Ts);

f2ss = frf2ss(G_uz2stage_frf, omegas, Nd2, ss_opts); % 12
sys_stage = f2ss.realize(ns2); % 12

figure(50); clf
pzplot(sys_stage);

Z = zero(sys_stage);
Z_eject = zpk([], Z(find(abs(Z) > 1)), 1, Ts);
Z_eject = Z_eject/dcgain(Z_eject);
sys_stage = minreal(Z_eject*sys_stage);

frfBode(sys_stage, freqs, F3, '--k', 'Hz');
plotPZ_freqs(sys_stage, F3);

% Now, we refine with the logrithmic least squares fit.
LGopts = optimoptions(@lsqnonlin, 'Display', 'iter',...
    'FunctionTolerance', 1e-7, 'MaxIter', 5000,'MaxFunctionEvaluations', 5000,...
    'StepTolerance', 1e-8, 'Jacobian','on', 'CheckGradients', false);

sos_fos = SosFos(sys_stage, 'iodelay', sys_stage.InputDelay);
LG = LogCostZPK(G_uz2stage_frf(1:k_estmax), omegas(1:k_estmax), sos_fos);
LG.solve_lsq(2, LGopts)
[sys_stage_log, p] = LG.sos_fos.realize();
sys_stage_log.InputDelay = max(round(p, 0), 0);
fprintf('LG says delay = %.2f\n', p);

frfBode(sys_stage_log, freqs, F4, '--k', 'Hz');
% plotPZ_freqs(sys_stage_log, F4);
% return
% Z = sort(zero(sys_stage_log))
% zej = Z(real(Z)<0.51)
% Z_eject = zpk([], zej, 1, Ts);
% Z_eject = Z_eject/dcgain(Z_eject);
% sys_stage_log = minreal(Z_eject*sys_stage_log);
% 
% sos_fos = SosFos(sys_stage_log, 'iodelay', p);
% LG = LogCostZPK(G_uz2stage_frf(1:k_estmax), omegas(1:k_estmax), sos_fos);
% LG.solve_lsq(2, LGopts)
% [sys_stage_log, p] = LG.sos_fos.realize();
% fprintf('LG says delay = %.2f\n', p);

sys_stage_log.InputDelay = max(round(p, 0), 0)+1;
frfBode(sys_stage_log, freqs, F4, '--b', 'Hz');

plotPZ_freqs(sys_stage_log, F4);
subplot(2,1,1)
title('Model with logfit')
figure(50); hold on;
pzplot(sys_stage_log, 'r')

%%
% ----------------------------------------------------------------
% --------------------- Now, Fit the drift model -----------------
addpath('hysteresis')
load(fullfile(PATHS.sysid, 'hysteresis', 'driftID_data_4-30-2018_01.mat'))
Ts = modelFit.frf.Ts;
G_uz2stage = sys_stage_log; %modelFit.models.G_uz2stage;



k_start = 2677;
y_exp = driftData.y_exp(k_start:end) - mean(driftData.y_exp(1:k_start));
u_exp = driftData.u_exp(k_start:end);
t_exp = (0:length(y_exp)-1)'*Ts;

theta0 = [0.9922    0.9997    0.9997    0.9927    .8];
lb = [-1, -1, -1, -1, -Inf];
ub = -lb;
np = 2;
normalize_dc = false;
Gvib = eject_gdrift(G_uz2stage, normalize_dc);

gdrift_cost = @(theta)fit_gdrift(theta, Gvib, y_exp, u_exp, t_exp, np);

opts = optimoptions(@lsqnonlin);
opts.MaxFunctionEvaluations = 5000;
opts.MaxIterations = 1000;
opts.Display = 'iter';
theta = lsqnonlin(gdrift_cost, theta0, lb, ub, opts);

gdrift = zpk(theta(np+1:end-1), theta(1:np), theta(end), Ts);

ydrift_est0 = lsim(Gvib*gdrift, u_exp, t_exp);
y_vib = lsim(Gvib, u_exp, t_exp);

figure(100)
clf;
h1 = plot(t_exp, y_exp);
h1.DisplayName = 'Exp. Step Response';
hold on
h2 = plot(t_exp, ydrift_est0, '-r');
h2.DisplayName = '$G_{vib}G_d$';


h3 = plot(t_exp, y_vib, ':k');
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

% saveas(F100, fullfile(PATHS.jfig, 'drift_fit.svg'));

%
% Now, see if we can get a bit better fit in the time domain by tweaking the
% first two complex resonances at 400 hz.
% figure(20); clf
% pzplot(Gvib);
% hold on
% 
% p = sort(pole(Gvib))
% z = sort(zero(Gvib))
% 
% p_eject = p(end-4+1:end)';
% z_eject = z(end-4+1:end)';
% 
% geject = zpk(z_eject, p_eject, 1, Ts);
% geject = geject/dcgain(geject);
% 
% Gvib_high = minreal(Gvib/geject)
% 
% pzplot(Gvib_high, 'r')
% 
% 
% 
% k0 = 1.0489;
% theta_0 = [real(z_eject([1,3])), imag(z_eject([1,3])),...
%            real(p_eject([1,3])), imag(p_eject([1,3])), k0]
% 
% g_fit = @(theta) zpk( [theta(1:2) + 1j*(theta(3:4)), theta(1:2) - 1j*(theta(3:4))],...
%           [theta(5:6)  + 1j*theta(7:8), theta(5:6) - 1j*theta(7:8)], theta(end), Ts);
% 
% 
% 
% yerr = @(theta) lsim(Gvib_high*g_fit(theta)*gdrift, u_exp, t_exp) - y_exp;
% 
% 
% theta = lsqnonlin(yerr, theta_0, lb, ub, opts);
% 
% 
% ydrift_est2 = lsim(Gvib_high*g_fit(theta)*gdrift, u_exp, t_exp);
 
% figure(100)
% plot(t_exp, ydrift_est2, 'g')
frfBode(Gvib_high*gdrift, freqs, F4, '--k', 'Hz')
%%
% ------------------- Fit Hysteresis -------------------------------------
hyst_file = 'hystID_data_5-4-2018_01.mat';
hyst_path = fullfile(PATHS.sysid, 'hysteresis', hyst_file);
load(hyst_path)

ux = hystData.u_exp;
yx = hystData.y_exp - mean(hystData.y_exp(1:500));
tvec = hystData.t_exp;

figure
plot(ux)

% start with prior infomation for guess
guess_data = load('steps_hyst_model.mat');

r = linspace(0, umax*0.5, nw);
r = ([0:nw-1]'./(nw) )*umax
d = [0; 4.5; 9];

d = [0; 5; 10]
w = ones(nw, 1);

umax = max(abs(ux));
nw = 7;

ws = [1; 0.0001; 0.0001];
G = Gvib*gdrift;
u_hyst_fun = @(theta) PIHyst.hyst_play_sat_op(ux, r, theta(1:nw), d, theta(nw+1:end), w*0);

cost_hyst_only = @(theta) downsample(lsim(G, u_hyst_fun(theta),... 
                tvec) - yx, 50); 
              
%%              
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
modelFit.models.G_uz2stage = sys_stage_log;
modelFit.models.G_uz2powI = G_deluz2Ipow*g_der;
modelFit.models.G_deluz2powI = G_deluz2Ipow;
modelFit.models.g_deluz2pow_1norm = nm1;
modelFit.models.du_max_nm1 = delumax;
modelFit.models.Gvib = Gvib;
modelFit.models.gdrift = gdrift;

if 0
    save(modelFit_file, 'modelFit');
end

