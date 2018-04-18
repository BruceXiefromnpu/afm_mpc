clc, clear
C_stage = 3.8e-6;
C1 = 4e-6;

Imax = 100e-3;
Ts = 40e-6;
MF_I_stage = load('pow_amp/FRF_data_current_stage2.mat');
MF_I_stage = MF_I_stage.modelFit;
G_uz2powI = MF_I_stage.models.G_uz2current1;
ms = 1e3;


R2 = 1.7e6;
R1 = 29.7e6;
Vdiv = R2/(R1 + R2);
Vdiv_gain = 1/Vdiv % from resistor measurement

load( fullfile(PATHS.exp, 'x-axis_sines_info_out_2-8-2018-01.mat'))
Gpow_hat = modelFit.models.G_uz2pow;

Gpow = ss(Gpow_hat*Vdiv_gain);
Gpow.InputDelay = 0;
Gpow = minreal(Gpow*zpk([], [0 0], 1, Ts))
Gpow = Gpow*5.8514/dcgain(Gpow)


dVmax = (Ts/C_stage)*Imax
dVmax = 0.5;

G_stage = modelFit.models.G_pow2uz/Vdiv_gain;

G_stage.InputDelay = 0;
sys_old = modelFit.models.G_uz2stage;
sys = G_stage*Gpow;
% sys = modelFit.models.G_uz2stage;
F1 = figure(1000);

[~,~,omegas] = bode(sys);
frfBode(Gpow, omegas/2/pi, F1,'-b', 'Hz');
frfBode(G_stage, omegas/2/pi, F1,'-k', 'Hz');
frfBode(sys_old, omegas/2/pi, F1,'-r', 'Hz');
frfBode(sys, omegas/2/pi, F1,'-g', 'Hz');

plotPZ_freqs(sys, F1);
plotPZ_freqs(G_stage, F1);

% figure(100)
% pzplot(Gpow, 'b', G_stage, 'k')
% 
% figure(3)
% pzplot(sys_old, 'g', sys, 'r')
% sys_nodelay = sys;
sys_nodelay = modelFit.models.G_uz2stage;
sys2 = modelFit.models.G_uz2stage;
[~, gdrift2] = eject_gdrift(sys2);

sys_nodelay.InputDelay = 0;
% [sys_nodelay, gdrift] = eject_gdrift(sys_nodelay)

%%
% Now, try the state constraint with the condensed formulation.
x0_pow = [0;0];

ref = 5; %OA works, QP fails
Nx = SSTools.getNxNu(sys_nodelay);
x0 = -Nx*ref;


% ----------------------------------------------------------------------- %
% -------------------- Now, try the time-optimal ------------------------ %

clc
k0 =78

C_effective = 4.0497e-05;
du_max = Ts*Imax/C_effective;

CON3 = CondenCon(Gpow, x0_pow, k0);
dVmax_opt = 5

CON3.add_state_con('box', dVmax_opt*01);
CON3.add_input_con('box', du_max);
CON3.update_binq

[sys_TO1, gdrift] = eject_gdrift(sys_nodelay);

sys_TO = SSTools.deltaUkSys(sys_TO1);
Nx_TO = SSTools.getNxNu(sys_TO);
%
x0_TO = Nx_TO*0;
xf = Nx_TO*ref*1;
P = path;
addpath(genpath(fullfile(getMatPath, 'solvers/cvx')))
C=[];

for k=0:k0-1
    C = [sys_TO.a^k*sys_TO.b C];
end
e_k0 = zeros(1,k0);
e_k0(end) = 1; % Kth unit vector. [0 0 .... 0 0 1]. To pick out last element of u.
I = eye(length(sys_TO.b));

cvx_begin
    variable u(k0);
    L = sys_TO.a^k0*x0_TO + C*u;
    minimize norm(L - xf)
    subject to
    [CON3.Ainq; -CON3.Ainq]*u <= [CON3.ubAinq; -CON3.lbAinq];
    u <= CON3.ub;
    -u<= -CON3.lb;
cvx_end
path(P);

fprintf('Results optval = %f\n', cvx_optval)
results.U = u;
results.cvx_optval = cvx_optval;
results

u = [u; repmat(0, 400, 1)];

t = [0:1:length(u)-1]'*Ts;
utraj_inv = lsim(gdrift, u, t);
utraj_inv = u;
uref = cumsum(utraj_inv);

sys_delay = sys_nodelay;
% y = lsim(sys_nodelay, uref, t, sys_nodelay.b*0);
y = lsim(sys_TO1, uref, t, sys_TO1.b*0);
sys_delay.InputDelay = 10;
ydelay = lsim(sys_delay, uref, t, sys_delay.b*0);
ypow = lsim(Gpow, utraj_inv, t, x0_pow*0);
Ipow = lsim(G_uz2powI, uref, t);

%
F1 = figure(9); clf;
subplot(3,1,1), hold on
h5 = plot(t, y, '-k');
h5.DisplayName = 'Time-Optimal';
grid on

% legend([h3,h4,h5]);
subplot(3,1,2), hold on
plot(t, cumsum(u), '-k')
title('u(k)')
grid on

subplot(3,1,3), hold on, grid on
plot(t, u, '-k');
plot(t(k0), u(k0), 'x');
title('$\Delta u_z$ (low voltage control)')
grid on

figure(200), clf
hold on
h6 = plot(t, ypow*C1/Ts);
h6.DisplayName = 'TO sim';
h6 = plot(t, Ipow, '-m');
% legend([h1, h2, h6])
% xlim([0, ])
    grid on
    title('Current, $I$')
    xlm = xlim;
    hold on;
    plot(xlm, [Imax, Imax], '--k')
    plot(xlm, -[Imax, Imax], '--k')



%%
% ---- Paths for shuffling data to labview and back. ------

% creat and pack data. Then save it. 
if 1
    slewfname_in = 'slewexp_datain_4-6-2018_01.csv';
    slewfpath_in = fullfile(pwd, slewfname_in)
    % slewfpath_in = fullfile(PATHS.exp, 'sysID', slewfname_in);

    slewfname_out = 'slewexp_dataout_4-6-2018_01.csv';
    % slewfpath_out = fullfile(PATHS.exp, 'sysID', slewfname_out);
    slewfpath_out = fullfile(pwd, slewfname_out);
    csvwrite(slewfpath_in, uref);
    clear vi;
    vipath = 'C:\Users\arnold\Documents\MATLAB\afm_mpc_journal\labview\play_nPoint_id_slew_OL_FIFO.vi'

    [e, vi] = setupVI(vipath, 'Abort', 0,...
                'umax', 9, 'data_out_path', slewfpath_out,...
                'traj_in_path', slewfpath_in, 'TsTicks', 1600);
    vi.Run
end


%%
% Now, read in data, and save to structure, and plot.

clc
data = csvread(slewfpath_out);
data = data(2:end, :);
size(data)

R_sense = 0.1;
ms = 1;
u_exp = data(:,1);
y_exp = data(:,2)- data(1,2);
upow_exp = data(:,3)*Vdiv_gain;
pow_I = data(:,4); 
pow_I = pow_I - mean(pow_I);
pow_I = pow_I/R_sense;
% amp output is centered around like 65 volts. Subtract that off so we can
% compare.
upow_exp = upow_exp - upow_exp(1); 

% upow_exp = data(:,3)*Vdiv_gain/dcgain(Gpow);
texp = [0:1:length(u_exp)-1]'*Ts;


figure(9);
subplot(3,1,1), hold on;
plot(texp*ms, y_exp, '--r')

% subplot(3,1,2)
% % plot(texp*ms, u_exp, '--r')
% 
% subplot(3,1,3), hold on

% plot(texp(1:end-1)*ms, du_exp*C1/Ts, '--k')
% title('(HV) $\Delta u(k)$', 'interpreter', 'latex')
% % ylim(1.2*[min(du_exp), max(du_exp)])

du_exp = diff(upow_exp);

figure(200)
plot(texp(1:end-1)*ms, du_exp*C1/Ts, '--r');
plot(texp*ms, pow_I, '--g')
title('current estimate')
grid on

%%
size(y_exp.Time)
figure(10), clf
subplot(2,1,1)
hold on
plot(t, ydelay)
plot(y_exp.Time, y_exp.Data, '--r')
grid on

subplot(2,1,2)
hold on
plot(t, uref)
plot(u_exp.Time, u_exp.Data, '--r')
grid on
legend('sim', 'exp')

figure(200)
dypow = diff(ypow_exp.Data)*Vdiv_gain;
h7 = plot(dypow(3:end), '--r');
h7.DisplayName = 'exp';

legend([h6, h7])

figure(400); clf
plot(diff(ypow))
hold on
plot(diff(dypow(3:end)))
title('$\Delta^2 y_{pow}$')
%%
sysinv = 1/zpk(sys_nodelay);
sysinv.InputDelay = 3;
sysinv = absorbDelay(sysinv)
ye = lsim(sysinv, y, t);
figure, plot(t, ye, t, y, '--')

%%


expOpts = stepExpOpts(linOpts, 'pstyle', 'k', 'name',  'AFM Stage');

afm_exp = stepExp(y_exp, u_exp, expOpts);
H2 = plot(afm_exp, F1);



