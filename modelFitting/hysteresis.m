clear
clc

modelFit_file = 'C:\Users\arnold\Documents\labview\sysID\data\x-axis_sines_info_intsamps_amp_1p0_out_4-2-2018offs_0-01.mat';
load(modelFit_file)

Ts = modelFit.frf.Ts;
G_uz2pow_rescaled = modelFit.models.G_uz2pow_rescaled;

R2 = 1.7e6;
R1 = 29.7e6;
Vdiv = R2/(R1 + R2);
Vdiv_gain = 1/Vdiv % from resistor measurement
% Vdiv_gain = 8.196/0.451; % measured from labview + 9v battery

C1 = 4e-6;
Imax = 100e-3;
dVmax_high = (Ts/C1)*Imax
dVmax_low = dVmax_high/dcgain(G_uz2pow_rescaled)

MF_I_stage = load('FRF_data_current_stage2.mat');
MF_I_stage = MF_I_stage.modelFit;
G_uz2powI = MF_I_stage.models.G_uz2current1;
ms = 1e3;

%% Case 1:
% The signal is not too aggressive. We expect nearly perfect match between
% power amp model and output. 
clc
To_half = 28*Ts;
To = 2*To_half;
wo = 2*pi/To;

t0 = [0:1:(3*To/Ts)-1]'*Ts;
du = 0.05;
du_vec = sign(sin(wo*t0))*du;
u_vec = [cumsum(du_vec); zeros(800, 1)];
t = (0:length(u_vec)-1)'*Ts;

F4 = figure(4); clf
subplot(3,1,1)
plot(t*ms, u_vec)
grid on
title('(LV) u(k)', 'interpreter', 'latex')
xlabel('time [ms]')

subplot(3,1,2)
plot(t(1:end-1)*ms, diff(u_vec))
grid on
title('(LV) $\Delta u(k)$', 'interpreter', 'latex')
xlabel('time [ms]')

% du_pow = lsim(Gpow/dcgain(Gpow), diff(u), t(1:end-1));
du_pow = lsim(G_uz2pow_rescaled, diff(u_vec), t(1:end-1));
Ipow = lsim(G_uz2powI, u_vec, t);


subplot(3,1,3)

plot(t(1:end-1)*ms, du_pow*C1/Ts)
grid on, hold on;
plot(t*ms, Ipow, '--r')
ylabel('"current"')
xlabel('time [ms]')
xlm = xlim;
plot(xlm, [dVmax_high, dVmax_high]*C1/Ts, ':k')
plot(xlm, -[dVmax_high, dVmax_high]*C1/Ts, ':k')
% ylim(1.2*[-dVmax_high, dVmax_high]*C1/Ts)

figure(5); clf
plot(t(1:end-1)*ms, du_pow*C1/Ts)
grid on, hold on;
plot(t*ms, Ipow, '--r')
ylabel('"current"')
xlabel('time [ms]')
xlm = xlim;
plot(xlm, [dVmax_high, dVmax_high]*C1/Ts, ':k')
plot(xlm, -[dVmax_high, dVmax_high]*C1/Ts, ':k')
% ylim(1.2*[-dVmax_high, dVmax_high]*C1/Ts)



slewfname_in = 'slewexp_datain_4-6-2018_01.csv';
slewfpath_in = fullfile(pwd, slewfname_in)
% slewfpath_in = fullfile(PATHS.exp, 'sysID', slewfname_in);

slewfname_out = 'slewexp_dataout_4-6-2018_01.csv';
% slewfpath_out = fullfile(PATHS.exp, 'sysID', slewfname_out);
slewfpath_out = fullfile(pwd, slewfname_out);
%%
% -----------------------RUN THE Experiment--------------------------------
if 1
    csvwrite(slewfpath_in, u_vec);
    clear vi;
    vipath = 'C:\Users\arnold\Documents\MATLAB\afm_mpc_journal\labview\play_nPoint_id_slew_OL_FIFO.vi'

    [e, vi] = setupVI(vipath, 'Abort', 0,...
                'umax', 9, 'data_out_path', slewfpath_out,...
                'traj_in_path', slewfpath_in, 'TsTicks', 1600);
    vi.Run
end
%%
data = csvread(slewfpath_out);
size(data)

R_sense = 0.1;

u_exp = data(:,1);
y_exp = data(:,2);
upow_exp = data(:,3)*Vdiv_gain;
pow_I = data(:,4); 
pow_I = pow_I - mean(pow_I);
pow_I = pow_I/R_sense;
% amp output is centered around like 65 volts. Subtract that off so we can
% compare.
upow_exp = upow_exp - upow_exp(1); 

% upow_exp = data(:,3)*Vdiv_gain/dcgain(Gpow);
texp = [0:1:length(u_exp)-1]'*Ts;

figure(F4);
subplot(3,1,1), hold on;
plot(texp*ms, u_exp, '--k')


subplot(3,1,3), hold on
du_exp = diff(upow_exp);
plot(texp(1:end-1)*ms, du_exp*C1/Ts, '--k')
title('(HV) $\Delta u(k)$', 'interpreter', 'latex')
% ylim(1.2*[min(du_exp), max(du_exp)])

subplot(3,1,3)
plot(texp*ms, pow_I, '--g')
title('current estimate')
grid on


figure(5)
plot(texp*ms, pow_I, '--g')
title('current estimate')
grid on
%%

%% Case 2:
% We try an aggressive, rapidly changing signal designed to exceed our
% bounds.

To_half = 20*Ts;
To = 2*To_half;
wo = 2*pi/To;

n0 = [0:To_half/Ts-1]';
t0 = n0*Ts;

% n1 = [n0(end)+1:n0(end)+1+3*(To/Ts)]';
n1 = [0:3.5*(To/Ts)]';
t1 = n1*Ts
du = 0.4;
du0 = sign(sin(t0*wo))*du;
du1 = sign(sin(0.5*wo*t1 + pi))*du;

du_vec = [du0; du1];

u_vec = [cumsum(du_vec); zeros(800, 1)];
t = (0:length(u_vec)-1)'*Ts;

F4 = figure(40); clf
subplot(3,1,1)
plot(t*ms, u_vec)
grid on
title('(LV) u(k)', 'interpreter', 'latex')
xlabel('time [ms]')

subplot(3,1,2)
plot(t(1:end-1)*ms, diff(u_vec))
grid on
title('(LV) $\Delta u(k)$', 'interpreter', 'latex')
xlabel('time [ms]')

% du_pow = lsim(Gpow/dcgain(Gpow), diff(u), t(1:end-1));
du_pow = lsim(G_uz2pow_rescaled, diff(u_vec), t(1:end-1));
Ipow = lsim(G_uz2powI, u_vec, t);

subplot(3,1,3)

plot(t(1:end-1)*ms, du_pow*C1/Ts)
grid on, hold on;
plot(t*ms, Ipow, '--r')
ylabel('"current"')
xlabel('time [ms]')
xlm = xlim;
plot(xlm, [dVmax_high, dVmax_high]*C1/Ts, ':k')
plot(xlm, -[dVmax_high, dVmax_high]*C1/Ts, ':k')
% ylim(1.2*[-dVmax_high, dVmax_high]*C1/Ts)

figure(50); clf
plot(t(1:end-1)*ms, du_pow*C1/Ts)
grid on, hold on;
plot(t*ms, Ipow, '--r')
ylabel('"current"')
xlabel('time [ms]')
xlm = xlim;
plot(xlm, [dVmax_high, dVmax_high]*C1/Ts, ':k')
plot(xlm, -[dVmax_high, dVmax_high]*C1/Ts, ':k')



slewfname_in = sprintf('slewexp_datain_4-2-2018_%0.2f.csv', max(diff(u_vec)));
slewfpath_in = fullfile(pwd, slewfname_in);
% slewfpath_in = fullfile(PATHS.exp, 'sysID', slewfname_in);

slewfname_out = sprintf('slewexp_dataout_4-2-2018_%0.2f.csv', max(diff(u_vec)));
% slewfpath_out = fullfile(PATHS.exp, 'sysID', slewfname_out);
slewfpath_out = fullfile(pwd, slewfname_out);
%%
% -----------------------RUN THE Experiment--------------------------------
if 1
    csvwrite(slewfpath_in, u_vec);
    clear vi;
    vipath = 'C:\Users\arnold\Documents\MATLAB\afm_mpc_journal\labview\play_nPoint_id_slew_OL_FIFO.vi'

    [e, vi] = setupVI(vipath, 'Abort', 0,...
                'umax', 9, 'data_out_path', slewfpath_out,...
                'traj_in_path', slewfpath_in, 'TsTicks', 1600);
    vi.Run
end
%%
data = csvread(slewfpath_out);
size(data)

u_exp = data(:,1);
y_exp = data(:,2);
upow_exp = data(:,3)*Vdiv_gain;
pow_I = data(:,4)/R_sense; 
pow_I = (pow_I - mean(pow_I)); %*1000;
% amp output is centered around like 65 volts. Subtract that off so we can
% compare.
upow_exp = upow_exp - upow_exp(1); 

% upow_exp = data(:,3)*Vdiv_gain/dcgain(Gpow);
texp = [0:1:length(u_exp)-1]'*Ts;
%
figure(40);
subplot(4,1,1), hold on;
plot(texp*ms, u_exp, '--k')


subplot(3,1,3), hold on
du_pow_exp = diff(upow_exp);
plot(texp(1:end-1)*ms, du_pow_exp*C1/Ts, '--')
title('(HV) $\Delta u(k)$', 'interpreter', 'latex')
ylim([1.2*min(du_pow_exp), 1.2*max(du_pow_exp)]*C1/Ts)


figure(50)
plot(texp*ms, pow_I, '--g')
title('current estimate')
grid on
%% Case 3:
% Try flattening the tops.
C1 = 4e-6;
Imax = 100e-3;
dVmax_high = ((1/C1)*Imax)/(1e6); % volts per micro-sec


To_half = 25*Ts;
To = 2*To_half;
wo = 2*pi/To;

n0 = [0:To_half/Ts-1]';
t0 = n0*Ts;

% n1 = [n0(end)+1:n0(end)+1+3*(To/Ts)]';
n1 = [0:3.5*(To/Ts)]';
t1 = n1*Ts;
du = 0.3;
du0 = sign(sin(t0*wo))*du;
du1 = sign(sin(0.5*wo*t1 + pi))*du;

du_vec = [du0; du1];

kcross = crossing(du_vec);
kcross = [kcross, length(du_vec)];
figure(4); clf
subplot(4,1,1)
hold on
plot(t(kcross)*ms, kcross*0, 'x')
duzero = zeros(40, 1);

du_vec_tmp = du_vec;
du_vec = [];
for k=1:length(kcross)-1
   du_vec = [du_vec; du_vec_tmp(kcross(k)+1:kcross(k+1)); duzero] ;
    
end

u_vec = [cumsum(du_vec); zeros(200, 1)];
t = (0:length(u_vec)-1)'*Ts;
%%
TO_dat = load('TO_fail.mat')
t = TO_dat.t;
u_vec = TO_dat.u;
du_vec = TO_dat.du;
R_sense = 0.1;

%%
F4 = figure(40); clf
subplot(4,1,1)
ax1 = gca;
plot(t*ms, u_vec)
grid on
title('(LV) u(k)', 'interpreter', 'latex')
xlabel('time [ms]')

subplot(4,1,2)
ax2 = gca;
plot(t(1:end-1)*ms, diff(u_vec))
grid on
title('(LV) $\Delta u(k)$', 'interpreter', 'latex')
xlabel('time [ms]')

% du_pow = lsim(Gpow/dcgain(Gpow), diff(u), t(1:end-1));
u_pow = lsim(G_uz2pow_rescaled, u_vec, t);
du_pow =1e-6* [diff(u_pow); 0]/Ts;
subplot(4,1,3)
ax3 = gca;
plot(t*ms, du_pow)
grid on, hold on;
xlabel('time [ms]')
xlm = xlim;
plot(xlm, [dVmax_high, dVmax_high], ':k')
plot(xlm, -[dVmax_high, dVmax_high], ':k')
% ylim(1.2*[-dVmax_high, dVmax_high])
ylim([1.2*min(du_pow), 1.2*max(du_pow)])
ylabel('du, [volts/$\mu$s]')


subplot(4,1,4)
ax4 = gca;
plot(t*ms, u_pow)
grid on;
xlabel('time [ms]')
ylabel('$y_{pow}$')
hold on

linkaxes([ax1, ax2, ax3, ax4], 'x')
xlim([0, t(end)*ms])


%
slewfname_in = sprintf('slewexp_datain_4-2-2018_%0.2f.csv', max(diff(u_vec)));
slewfpath_in = fullfile(pwd, slewfname_in);
% slewfpath_in = fullfile(PATHS.exp, 'sysID', slewfname_in);

slewfname_out = sprintf('slewexp_dataout_4-2-2018_%0.2f.csv', max(diff(u_vec)));
% slewfpath_out = fullfile(PATHS.exp, 'sysID', slewfname_out);
slewfpath_out = fullfile(pwd, slewfname_out);
%
% -----------------------RUN THE Experiment--------------------------------
if 1
    csvwrite(slewfpath_in, u_vec);
    clear vi;
    vipath = 'C:\Users\arnold\Documents\MATLAB\afm_mpc_journal\labview\play_nPoint_id_slew_OL_FIFO.vi';

    [e, vi] = setupVI(vipath, 'Abort', 0,...
                'umax', 9, 'data_out_path', slewfpath_out,...
                'traj_in_path', slewfpath_in, 'TsTicks', 1600);
    vi.Run
end
%
data = csvread(slewfpath_out);
size(data)

u_exp = data(:,1);
y_exp = data(:,2);
upow_exp = data(:,3)*Vdiv_gain;
% amp output is centered around like 65 volts. Subtract that off so we can
% compare.
upow_exp = upow_exp - upow_exp(1); % - mean(upow_exp(end-10:end)); 

pow_I = data(:,4)/R_sense; %/VI_sense_gain/R_sense;
% pow_I = (pow_I - mean(pow_I)); %*1000;

% upow_exp = data(:,3)*Vdiv_gain/dcgain(Gpow);
texp = (0:1:length(u_exp)-1)'*Ts;
%
figure(F4);
subplot(4,1,1), hold on;
plot(texp*ms, u_exp, '--k')


subplot(4,1,3), hold on
du_pow_exp =1e-6* diff(upow_exp)/Ts;
plot(texp(1:end-1)*ms, du_pow_exp, '--')
title('(HV) $\Delta u(k)$', 'interpreter', 'latex')
ylim([1.2*min(du_pow_exp), 1.2*max(du_pow_exp)])

subplot(4,1,4)
plot(texp(1:end-1)*ms, upow_exp(2:end), '--')
title('(HV) $y_{pow}$', 'interpreter', 'latex')
grid on

figure(50)
% subplot(5,1,5)
plot(texp*ms, pow_I)
title('current estimate')
grid on
ax5 = gca
linkaxes([ax1, ax2, ax3, ax4, ax5], 'x')




