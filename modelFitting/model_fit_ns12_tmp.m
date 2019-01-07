clear
clc

saveon = false;

addpath(fullfile(PATHS.step_exp(), 'functions'))
addpath(fullfile(PATHS.step_exp(), 'functions', 'canon'))
if ~ispc
addpath('~/gradschool/sysID/matlab/functions')
else
  addpath('C:\Users\arnold\Documents\labview\sysID\matlab\functions\')
end
% modelFit_file = 'FRF_data_current_stage.mat';
% data_fname = 'x-axis_sines_infoFourierCoef_5-30-2018-01.mat';
data_fname = 'x-axis_sines_infoFourierCoef_9-10-2018-02.mat';
modelFit_file = fullfile(PATHS.sysid, 'FRF_data', data_fname);

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
frfBode(G_uz2stage_frf, freqs, F3, 'Hz', 'r');
subplot(2,1,1)
title('ERA fit')
frfBode(G_uz2stage_frf, freqs, F4,  'Hz', 'r');
subplot(2,1,1)
title('log fit')

Nd2 = 10;
ns2 = 12;
k_estmax = 273+75;
ss_opts = frf2ss_opts('Ts', Ts);

f2ss = frf2ss(G_uz2stage_frf, omegas, Nd2, ss_opts); % 12
sys_stage = f2ss.realize(ns2); % 12
%%%
plants = CanonPlants.plants_ns14();
sys_stage = plants.G_uz2stage;
%%%
%         
Z = zero(sys_stage);
Z_eject = zpk([], Z(find(abs(Z) > 1)), 1, Ts);
Z_eject = Z_eject/dcgain(Z_eject)
sys_stage = minreal(Z_eject*sys_stage);

frfBode(sys_stage, freqs, F3, 'Hz', '--k');
plotPZ_freqs(sys_stage, F3);


figure(F4);
subplot(2,1,1)
ylm = ylim;
plot([freqs(k_estmax), freqs(k_estmax)], ylm, 'k');

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

frfBode(sys_stage_log, freqs, F4,  'Hz', '--k');
plotPZ_freqs(sys_stage_log, F4);

[plants2, frf_data2] = CanonPlants.plants_ns14();
frfBode(frf_data2.G_uz2stage, frf_data2.freqs_Hz, F4, 'Hz', 'b')

% frfBode(plants2.G_uz2stage, frf_data2.freqs_Hz, F4, 'Hz', '--g')

figure(50);
pzplot(sys_stage_log);

% return
% Z = sort(zero(sys_stage_log));
% Z_eject = zpk([], Z(1:1), 1, Ts);
% Z_eject = Z_eject/dcgain(Z_eject);
% sys_stage_log = minreal(Z_eject*sys_stage_log);
% 
% sos_fos = SosFos(sys_stage_log, 'iodelay', p);
% LG = LogCostZPK(G_uz2stage_frf(1:k_estmax), omegas(1:k_estmax), sos_fos);
% LG.solve_lsq(2, LGopts)
% [sys_stage_log, p] = LG.sos_fos.realize();
% fprintf('LG says delay = %.2f\n', p);
% 
% sys_stage_log.InputDelay = max(round(p, 0), 0);
% frfBode(sys_stage_log, freqs, F4, 'Hz', '--b');
% 
% plotPZ_freqs(sys_stage_log, F4);
%
%%
% ---------------------------------------------------------------- %
% --------- Second, we work on the powI system -------------------- %

% Visualize everything
F10 = figure(10); clf
frfBode(G_uz2powI_frf, freqs, F10, 'Hz', 'r');


% Divide the derivative (z-1) out of the systems FRF.
ejw = exp(1j*Ts*modelFit.frf.w_s(:));

% der_frf_ct = (1j*modelFit.frf.w_s(:));
derz_frf = ejw - 1;
G_deluz2powI_frf = G_uz2powI_frf./derz_frf;

mag_max = max(abs(G_deluz2powI_frf));
derz_frf_adj = derz_frf*mag_max;
I_max = 0.1; %Amps
fprintf('(Hinf)Mag-max = %.3f, psudeo deltaUk_max = %.3f\n', mag_max, I_max/mag_max);

F20 = figure(20); clf;
frfBode(G_deluz2powI_frf, freqs, F20, 'Hz', 'r');
%
Nd1 = 4;
ss_opts = frf2ss_opts('Ts', Ts);

f2ss = frf2ss(G_deluz2powI_frf, omegas, Nd1, ss_opts); % 12
sys = f2ss.realize(8); 
% Remove NMP zeros
% Z = zero(sys);
% Z_eject = zpk([], Z(find(abs(Z) > 1)), 1, Ts);
% Z_eject = Z_eject/dcgain(Z_eject);
% sys = minreal(Z_eject*sys)


g_der = zpk([1], [], 1, Ts); % 
frfBode(sys, freqs, F20,  'Hz', '--k');
sys.InputDelay = 4;
frfBode(sys*g_der, freqs, F10, 'Hz', '--k');

sos_fos = SosFos(sys, 'iodelay', sys.InputDelay);
LG = LogCostZPK(G_deluz2powI_frf, freqs*2*pi, sos_fos);
LG.solve_lsq(2, LGopts)

[G_deluz2Ipow, p] = LG.sos_fos.realize();
G_deluz2Ipow.InputDelay = max(round(p, 0), 0);

frfBode(G_deluz2Ipow*g_der, freqs, F10, 'Hz', '--b');
frfBode(G_deluz2Ipow, freqs, F20,  'Hz', '--b');

plotPZ_freqs(G_deluz2Ipow*g_der, F10);

% compute impulse response for 1-norm:
[y, t] = impulse(G_deluz2Ipow);
y = y*Ts; % matlab makes the impulse 1/Ts tall.
nm1 = sum(abs(y))

delumax = I_max/nm1;
fprintf('(Hinf)Mag-max = %.3f, psudeo deltaUk_max = %.3f\n', mag_max, I_max/mag_max);

fprintf('(BIBO) ||G_delu2Ipow||_1 = %.3f, deltaUk_max = %.3f\n', nm1, delumax);


%
% ----------------------------------------------------------------
% --------------------- Now, Fit the drift model -----------------
addpath('hysteresis')
%%
load(fullfile(PATHS.sysid, 'hysteresis', 'driftID_data_5-29-2018_01.mat'))
% load(fullfile(PATHS.sysid, 'hysteresis', 'driftID_data_09-10-2018_01_amp_0p15.mat'))

% % % load(fullfile(PATHS.sysid, 'hysteresis', 'driftID_data_06-05-2018_01_amp_1p0.mat'))


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
normalize_dc = true;
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

figure(105)
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
% ylim([-0.005, 0.15])
grid on
% grid minor;
xlabel('time [s]')
ylabel('$y_X$ [v]')
ax = gca;
% set(ax, 'XTick', (0:0.05:0.3), 'YTick', (0:0.025:0.15))


%%
% ------------------- Fit Hysteresis + Sat -------------------------------------

fprintf('============================================\n')
% hyst_file = 'hystID_data_5-4-2018_01.mat';
hyst_file = 'hystID_data_10-Sep-2018_01.mat';
% % % hyst_file = 'hystID_data_27-Aug-2018_01.mat';
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
[rp, wp] = PIHyst.invert_hyst_PI(r, w);
[dp, wsp] = PIHyst.invert_sat(d, ws);
hyst_sat = struct('r', r, 'w', w, 'rp', rp, 'wp', wp,...
                  'd', d, 'ws', ws, 'dp', dp, 'wsp', wsp);

fprintf('r w d ws\n')                
[r, w, d, ws]

clear PIHyst
[r2, w2] = PIHyst.fit_hyst_weights(downsample(ux, 100), downsample(yprime, 100), Nhyst);

u_pisat = PIHyst.hyst_play_sat_op(ux, r, w, d, ws, w*0);
u_pi  = PIHyst.hyst_play_op(ux, r2, w2, w*0);

figure(1002); clf
plot(tvec, yprime)
hold on
plot(tvec, u_pisat)
grid on

y_hyst_sat = lsim(Gvib*gdrift, u_pisat, tvec);
y_hyst = lsim(Gvib*gdrift, u_pi, tvec);

figure(2002); clf
plot(tvec, yx)
hold on
plot(tvec, y_hyst_sat)
plot(tvec, y_hyst, '--')
legend('y-exp', 'y-hyst-sat', 'y-hyst')
grid on

[rp2, wp2] = PIHyst.invert_hyst_PI(r2, w2);

hyst = struct('r', r2, 'w', w2, 'rp', rp2, 'wp', wp2);
%%
frfBode(gdrift*Gvib, freqs, F4, 'Hz', '--m')
%%
plotPZ_freqs(gdrift*Gvib, F4);
%%
% modelFit.models.G_uz2powI = G_deluz2Ipow*g_der;
% modelFit.models.G_deluz2powI = G_deluz2Ipow;
% modelFit.models.g_deluz2pow_1norm = nm1;
% modelFit.models.du_max_nm1 = delumax;
modelFit.models.G_uz2stage = sys_stage_log;
modelFit.models.G_uz2powI = G_deluz2Ipow*g_der;
modelFit.models.G_deluz2powI = G_deluz2Ipow;
modelFit.models.g_deluz2pow_1norm = nm1;
modelFit.models.du_max_nm1 = delumax;
modelFit.models.Gvib = Gvib;
modelFit.models.gdrift = gdrift;
modelFit.models.hyst = hyst;
modelFit.models.hyst_sat = hyst_sat;
if saveon
    save(modelFit_file, 'modelFit');
end












