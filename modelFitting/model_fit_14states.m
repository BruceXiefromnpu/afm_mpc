clear
clc
% The goal of this script is to get a better fit at the first complex
% resonance. It seems we need 14 states to do that.

addpath('~/gradschool/sysID/matlab/functions')

% modelFit_file = fullfile(PATHS.sysid, 'FRF_data_current_stage2.mat');
modelFit_file = fullfile(PATHS.sysid,'FRF_data', 'x-axis_sines_info_HIRESFourierCoef_5-14-2018-01.mat');

load(modelFit_file)

freqs = modelFit.frf.freq_s(:);

% G_uz2powI_frf = modelFit.frf.G_uz2powI;
G_uz2stage_frf = modelFit.frf.G_uz2stage;
G_uz2powI_frf = modelFit.frf.G_uz2powI;

omegas = freqs*2*pi;
Ts = modelFit.frf.Ts;

%
% --------------------------------------------------------------- %
% --------- First, we work on the stage system ----------------- %
% 
F3 = figure(3); clf
F4 = figure(4); clf
frfBode(G_uz2stage_frf, freqs, F3, 'r', 'Hz');
frfBode(G_uz2stage_frf, freqs, F4, 'r', 'Hz');

Nd2 = 10;
ns2 = 14;
k_estmax = 250;
ss_opts = frf2ss_opts('Ts', Ts);

f2ss = frf2ss(G_uz2stage_frf, omegas, Nd2, ss_opts); % 12
sys_stage = f2ss.realize(ns2); % 12
                               
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
plotPZ_freqs(sys_stage_log, F4);

figure(50);
pzplot(sys_stage_log);

% return
%
Z = sort(zero(sys_stage_log));
Z_eject = zpk([], Z(abs(Z)>1), 1, Ts);
Z_eject = Z_eject/dcgain(Z_eject);
sys_stage_log = minreal(Z_eject*sys_stage_log);

sos_fos = SosFos(sys_stage_log, 'iodelay', p);
LG = LogCostZPK(G_uz2stage_frf(1:k_estmax), omegas(1:k_estmax), sos_fos);
LG.solve_lsq(2, LGopts)
[sys_stage_log, p] = LG.sos_fos.realize();
fprintf('LG says delay = %.2f\n', p);

sys_stage_log.InputDelay = max(round(p, 0), 0);
frfBode(sys_stage_log, freqs, F4, '--b', 'Hz');

plotPZ_freqs(sys_stage_log, F4);

figure(50); hold on;
pzplot(sys_stage_log, 'r')
%%

% ---------------------------------------------------------------- %
% --------- Second, we work on the powI system -------------------- %

% Visualize everything
F10 = figure(10); clf
frfBode(G_uz2powI_frf, freqs, F10, 'r', 'Hz');


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
frfBode(G_deluz2powI_frf, freqs, F20, 'r', 'Hz');
%
Nd1 = 4;
ss_opts = frf2ss_opts('Ts', Ts);

f2ss = frf2ss(G_deluz2powI_frf, omegas, Nd1, ss_opts); % 12
sys = f2ss.realize(8); 
% Remove NMP zeros
Z = zero(sys);
Z_eject = zpk([], Z(find(abs(Z) > 1)), 1, Ts);
Z_eject = Z_eject/dcgain(Z_eject);
sys = minreal(Z_eject*sys)


g_der = zpk([1], [], 1, Ts); % 
frfBode(sys, freqs, F20, '--k', 'Hz');
sys.InputDelay = 3;
frfBode(sys*g_der, freqs, F10, '--k', 'Hz');

sos_fos = SosFos(sys, 'iodelay', sys.InputDelay);
LG = LogCostZPK(G_deluz2powI_frf, freqs*2*pi, sos_fos);
LG.solve_lsq(2, LGopts)

[G_deluz2Ipow, p] = LG.sos_fos.realize();
G_deluz2Ipow.InputDelay = max(round(p, 0), 0);

frfBode(G_deluz2Ipow*g_der, freqs, F10, '--b', 'Hz');
frfBode(G_deluz2Ipow, freqs, F20, '--b', 'Hz');

plotPZ_freqs(G_deluz2Ipow*g_der, F10);

% compute impulse response for 1-norm:
[y, t] = impulse(G_deluz2Ipow);
y = y*Ts; % matlab makes the impulse 1/Ts tall.
nm1 = sum(abs(y))

delumax = I_max/nm1;
fprintf('(Hinf)Mag-max = %.3f, psudeo deltaUk_max = %.3f\n', mag_max, I_max/mag_max);

fprintf('(BIBO) ||G_delu2Ipow||_1 = %.3f, deltaUk_max = %.3f\n', nm1, delumax);


%%

modelFit.models.G_uz2stage = sys_stage_log;
modelFit.models.G_uz2powI = G_deluz2Ipow*g_der;
modelFit.models.G_deluz2powI = G_deluz2Ipow;
modelFit.models.g_deluz2pow_1norm = nm1;
modelFit.models.du_max_nm1 = delumax;
if 1
    save(modelFit_file, 'modelFit');
end

%%
% Now, Fit the drift model.
addpath('hysteresis')
load(fullfile(PATHS.sysid, 'hysteresis', 'driftID_data_4-30-2018_01.mat'))
Ts = modelFit.frf.Ts;
G_uz2stage = sys_stage_log; %modelFit.models.G_uz2stage;



k_start = 2677;
y_exp = driftData.y_exp(k_start:end) - mean(driftData.y_exp(1:k_start));
u_exp = driftData.u_exp(k_start:end);
t_vec = (0:length(y_exp)-1)'*Ts;

% F100 = figure(100); clf
% plot_ind(t_vec, y_exp);


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

% clf;
% h1 = plot(t_vec, y_exp);
% h1.DisplayName = 'Exp. Step Response';
figure(100)
hold on
h2 = plot(t_vec, ydrift_est0, '-g');
h2.DisplayName = '$G_{vib}G_d$';

%%
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
% saveas(F100, fullfile(PATHS.jfig, 'drift_fit.svg'));








