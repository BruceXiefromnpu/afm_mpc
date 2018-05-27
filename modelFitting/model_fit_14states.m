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
% frfBode(Gvib_high*gdrift, freqs, F4, '--k', 'Hz')
%%
% ------------------- Fit Hysteresis + Sat -------------------------------------


hyst_file = 'hystID_data_5-4-2018_01.mat';
hyst_path = fullfile(PATHS.sysid, 'hysteresis', hyst_file);
load(hyst_path)

kk = length(hystData.y_exp);
ux = hystData.u_exp(1:kk);
yx = hystData.y_exp(1:kk) - mean(hystData.y_exp(1:500));
tvec = hystData.t_exp(1:kk);

umax = max(abs(ux));
ymax = max(abs(yx));

figure(500); clf
plot(ux)
grid on

Nhyst = 7;
nw = Nhyst;
Nsat = 5;

yprime = lsim(1/(gdrift*dcgain(Gvib)), yx, tvec);
[r, w, d, ws] = PIHyst.fit_hyst_sat_weights(downsample(ux, 100), downsample(yprime, 100), Nhyst, Nsat);
[rp, wp] = PIHyst.invert_hyst_PI(r, w);
[dp, wsp] = PIHyst.invert_sat(d, ws);
hyst_sat = struct('r', r, 'w', w, 'rp', rp, 'wp', wp,...
                  'd', d, 'ws', ws, 'dp', dp, 'wsp', wsp);


%%
clear PIHyst
[r2, w2] = PIHyst.fit_hyst_weights(downsample(ux, 100), downsample(yprime, 100), Nhyst);

u_pisat = PIHyst.hyst_play_sat_op(ux, r, w, d, ws, w*0);
u_pi  = PIHyst.hyst_play_op(ux, r2, w2, w*0);

figure(1000); clf
plot(tvec, yprime)
hold on
plot(tvec, u_pisat)
grid on

y_hyst_sat = lsim(Gvib*gdrift, u_pisat, tvec);
y_hyst = lsim(Gvib*gdrift, u_pi, tvec);

figure(2000); clf
plot(tvec, yx)
hold on
plot(tvec, y_hyst_sat)
plot(tvec, y_hyst, '--')
grid on

[rp2, wp2] = PIHyst.invert_hyst_PI(r2, w2)

hyst = struct('r', r2, 'w', w2, 'rp', rp2, 'wp', wp2)

%%
modelFit.models.G_uz2stage = sys_stage_log;
modelFit.models.G_uz2powI = G_deluz2Ipow*g_der;
modelFit.models.G_deluz2powI = G_deluz2Ipow;
modelFit.models.g_deluz2pow_1norm = nm1;
modelFit.models.du_max_nm1 = delumax;
modelFit.models.Gvib = Gvib;
modelFit.models.gdrift = gdrift;
modelFit.models.hyst = hyst;
modelFit.models.hyst_sat = hyst_sat;
if 1
    save(modelFit_file, 'modelFit');
end


return
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

r = ([0:nw-1]'./(nw) )*umax;

%n_d = 0;
id_plus = (1:n_d);
id_neg = (-n_d:-1);
dplus = ((id_plus - 0.5)/n_d ) * ymax;
dmin = ( (id_neg + 0.5)/n_d ) *ymax;
dp = [dmin, 0, dplus]';

% create dummy vectors:
w = r*0;
ws = dp*0;

[~, HU_mat] = PIHyst.hyst_play_op(downsample(ux, 100), r, r*0, r*0);


[~, Syp_mat] = PIHyst.sat_op(downsample(yprime, 100), dp, dp*0);

H = [HU_mat'; -Syp_mat']*[HU_mat, -Syp_mat];

UH = -eye(size(w,1));
neg1 = -ones(1, n_d);
N = 2*n_d+1;
US = zeros(N, N);
On = ones(n_d+1);
Us1 = -triu(On);
Us2 = -triu(On)';
US = blkdiag(Us1(1:n_d, 1:n_d), Us2);
US(1:n_d, n_d+1) = -1;

Ainq = blkdiag(UH, US);
binq = [w(:)*0; ws(:)*0] - .011;

Aeq = [ymax*ones(1, length(r)) - r', zeros(1, length(dp))];
beq = ymax;

opts = optimset('quadprog');
opts.Display = 'off';
[w_wsp, JVAL] = quadprog(H, H(:,1)*0, Ainq, binq, Aeq, beq, [], [], [], opts);
fprintf('N-r = %.0f,  N-d = %.0f, JVAL = %g\n', length(r), length(d), JVAL);

w = w_wsp(1:length(r));
wsp = w_wsp(length(r)+1:end);

[d, ws] = PIHyst.invert_sat(dp, wsp);
