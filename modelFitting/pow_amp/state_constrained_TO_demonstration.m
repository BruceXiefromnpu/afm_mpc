clc, clear
C_stage = 3.8e-6;
Imax = 100e-3;
Ts = 40e-6;



% This loads up our models of the TF from uz to output current.
MF_GpowI_stage = load('pow_amp/FRF_data_current_stage2.mat')
MF_GpowI_cap = load('pow_amp/FRF_data_current_caponly2.mat')
G_uz2powI_stage = MF_GpowI_stage.modelFit.models.G_uz2current1;
G_uz2powI_cap = MF_GpowI_cap.modelFit.models.G_uz2current1;

%%

% This loads up the models we were using before. Those were un-scaled, so we 
% must scale them here.
load( fullfile(PATHS.exp, 'x-axis_sines_info_out_2-8-2018-01.mat'))
Vdiv = 2.56/(19.61+2.56);
% Vdiv_gain = 1/Vdiv; % from resistor measurement
Vdiv_gain = 10.6; % From 9v battery measurement

Gpow_hat = modelFit.models.G_uz2pow;

Gpow = ss(Gpow_hat*Vdiv_gain);
Gpow.InputDelay = 0;
Gpow = minreal(Gpow*zpk([], [0 0], 1, Ts))


dVmax = (Ts/C_stage)*Imax
dVmax = 0.5;

G_stage = modelFit.models.G_pow2uz/Vdiv_gain;

G_stage.InputDelay = 0;
sys_old = modelFit.models.G_uz2stage;
sys = G_stage*Gpow;
% sys = modelFit.models.G_uz2stage;
F1 = figure(1);

[~,~,omegas] = bode(sys);
frfBode(Gpow, omegas/2/pi, F1,'Hz', '-b');
frfBode(G_stage, omegas/2/pi, F1,'Hz', '-k');
frfBode(sys_old, omegas/2/pi, F1,'Hz', '-r');
frfBode(sys, omegas/2/pi, F1,'Hz', '-g');

plotPZ_freqs(sys, F1)
plotPZ_freqs(G_stage, F1)


%%
% sys_nodelay = sys;
sys_nodelay = modelFit.models.G_uz2stage;
sys_nodelay.InputDelay = 0;


% ----------------------------------------------------------------------- %
% -------------------- Now, try the time-optimal ------------------------ %
clc
k0 = 108

sys_recyc = SSTools.deltaUkSys(sys_nodelay);
[Nx, Nu] = SSTools.getNxNu(sys_recyc);
x0_pow = [0;0];
du_max = 2.0;
ref = 4; %OA works, QP fails
x0 = -Nx*ref;

CON3 = CondenCon(Gpow, x0_pow, k0);
CON3.add_state_con('box', dVmax);
CON3.add_input_con('box', du_max);
CON3.update_binq
% f_1 = CondensedTools.init_cond_resp_matrix(PHI_pow, 1, k0, C_pow);
% H_1 = CondensedTools.zero_state_output_resp(Gpow, k0, 1);
% 
% binq = ones(k0,1)*dVmax - f_1*x0_pow*0;
sys_TO = eject_gdrift(sys_nodelay);

sys_TO = SSTools.deltaUkSys(sys_TO);
Nx_TO = SSTools.getNxNu(sys_TO);

x0_TO = Nx_TO*0;
xf = Nx_TO*ref;

addpath(genpath(fullfile(getMatPath, 'solvers/cvx')))
C=[];
% du_max = 0.2;
for k=0:k0-1
    C = [sys_TO.a^k*sys_TO.b C];
end
e_k0 = zeros(1,k0);
e_k0(end) = 1; % Kth unit vector. [0 0 .... 0 0 1]. To pick out last element of u.
I = eye(length(sys_TO.b));

cvx_begin
    variable u(k0);
%     variable t;
    L = sys_TO.a^k0*x0_TO + C*u;
    minimize norm(L - xf)
    % minimize norm(u)
    % minimize norm(t)
    subject to
%     norm(u, Inf) <= du_max;
    [CON3.Ainq; -CON3.Ainq]*u <= [CON3.ubAinq; -CON3.lbAinq];
    u <= CON3.ub;
    -u<= -CON3.lb;
    % Enforce steady state requirement
    % xss = Axss + Buss --> (I-A)xss = Buss
    % (I-A)xss - Buss = 0. Not sure this is necessary.
    % If minimization succeeds, should be the case anyway.
    % (I-sys_TO.a)*L - sys_TO.b*(e_k0*u) ==0;
    % L - xf == 0;
cvx_end
rmpath(genpath(fullfile(getMatPath, 'solvers/cvx')))

fprintf('Results optval = %f\n', cvx_optval)
results.U = u;
results.cvx_optval = cvx_optval;
results

u = [u; repmat(0, 400, 1)];

t = [0:1:length(u)-1]'*Ts;
y = lsim(sys_TO, u, t, x0_TO*0);
ypow = lsim(Gpow, u, t, x0_pow*0);


I_current_stage = lsim(G_uz2powI_stage, u, t);
I_current_cap = lsim(G_uz2powI_cap, u, t);
%%
figure(9)
subplot(3,1,1), hold on
h5 = plot(t, y, '-.g');
h5.DisplayName = 'Time-Optimal'
legend([h5]);
subplot(3,1,2), hold on
plot(t, u, '-.g')
plot(t(k0), u(k0), 'x')
title('$\Delta u$')



subplot(3,1,3), hold on
plot(t, cumsum(u), '-.g')

title('u(k)')

figure(200); clf
hold on
h6 = plot(t, ypow*C_stage/Ts);
h6.DisplayName = '$\Delta V$ estimate current';

h7 = plot(t, I_current_stage);
h7.DisplayName = 'Simulated Current (model ID)';

h8 = plot(t, I_current_cap, '--k');
h8.DisplayName = 'Simulated current, pure capacitance load model';

legend([h6, h7, h8])
xlim([0,t( 300)])
xlm = xlim
plot(xlm, [Imax, Imax], ':k', 'LineWidth', 1.5)
plot(xlm, -[Imax, Imax], ':k', 'LineWidth', 1.5)
grid on

ylabel('amps')

xlim([0, 0.006])
%%
rmpath(genpath(fullfile(getMatPath, 'solvers/cvx')))






