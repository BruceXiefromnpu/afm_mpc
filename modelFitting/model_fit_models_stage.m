clear
clc

% modelFit_file = 'FRF_data_current_stage.mat';
modelFit_file = 'pow_amp/FRF_data_current_stage2.mat';

load(modelFit_file)

freqs = modelFit.frf.freq_s(:);

% G_uz2powI_frf = modelFit.frf.G_uz2powI;
G_uz2stage_frf = modelFit.frf.G_uz2stage;
G_uz2powI_frf = modelFit.frf.G_uz2powI;

omegas = freqs*2*pi;
Ts = modelFit.frf.Ts;

%%
% --------------------------------------------------------------- %
% --------- Second, we work on the stage system ----------------- %
% 
F3 = figure(3); clf
F4 = figure(4); clf
frfBode(G_uz2stage_frf, freqs, F3, 'r', 'Hz');
frfBode(G_uz2stage_frf, freqs, F4, 'r', 'Hz');

Nd2 = 10;
ns2 = 12;
k_estmax = 215;
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

frfBode(sys_stage_log, freqs, F4, '--k', 'Hz')
plotPZ_freqs(sys_stage_log, F4)

figure(50);
pzplot(sys_stage_log);

return

Z = sort(zero(sys_stage_log));
Z_eject = zpk([], Z(1:1), 1, Ts);
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
%%

% Visualize everything
F1 = figure(1); clf
frfBode(G_uz2powI_frf, freqs, F1, 'r', 'Hz');




% ---------------------------------------------------------------- %
% --------- First, we work on the powI system -------------------- %
% Divide the derivative (z-1) out of the systems FRF.
ejw = exp(1j*Ts*modelFit.frf.w_s(:));

der_frf_ct = (1j*modelFit.frf.w_s(:));
der_frf = ejw - 1;
G_uz2powI_int_frf = G_uz2powI_frf./der_frf;


F2 = figure(2); clf;
frfBode(G_uz2powI_int_frf, freqs, F2, 'r', 'Hz');

Nd1 = 4;
ss_opts = frf2ss_opts('Ts', Ts);

f2ss = frf2ss(G_uz2powI_int_frf, omegas, Nd1, ss_opts); % 12
sys = f2ss.realize(12); % 12

g_der = zpk([1], [], 1, Ts); % 
frfBode(sys, freqs, F2, '--k', 'Hz');
sys.InputDelay = 3;
frfBode(sys*g_der, freqs, F1, '--k', 'Hz');
%

sos_fos = SosFos(sys, 'iodelay', sys.InputDelay);
LG = LogCostZPK(G_uz2powI_int_frf, freqs*2*pi, sos_fos);
LG.solve_lsq(2, LGopts)

[sys3, p] = LG.sos_fos.realize();
sys3.InputDelay = max(round(p, 0), 0);
frfBode(sys3*g_der, freqs, F1, '--b', 'Hz');
frfBode(sys3, freqs, F2, '--b', 'Hz');

plotPZ_freqs(sys3*g_der, F1);





%%


modelFit.models.G_uz2powI1 = sys3*g_der;
if 0
    save(modelFit_file, 'modelFit');
end












