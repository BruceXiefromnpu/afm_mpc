clear
clc

% modelFit_file = '/media/labserver/mpc-journal/x-axis_sines_info_out_2-8-2018-01.mat';
modelFit_file = 'C:\Users\arnold\Documents\labview\sysID\data\x-axis_sines_info_intsamps_amp_1p0_out_4-2-2018offs_0-01.mat';

load(modelFit_file)
% load /media/labserver/mpc-journal/sysID/x-axis_sines_info_matchTs_out_2-14-2018-01.mat
whos

S_ij = modelFit.frf.S_ij;
freqs = modelFit.frf.freq_s;



P_pow2stage = squeeze(S_ij(1,2, :)./S_ij(2,2, :));
P_uz2stage = modelFit.frf.G_uz2stage;
P_uz2pow = modelFit.frf.G_uz2pow;
%%

% Visualize everything
F1 = figure(1); clf
frfBode(P_uz2pow, freqs, F1, 'Hz', 'b');

omegas = freqs*2*pi;
Ts = modelFit.frf.Ts;
P_uz2pow_frd = frd(P_uz2pow, omegas, Ts);

G_uz2pow = tfest(P_uz2pow_frd, 2, 0,5, 'Ts', Ts);
% G_uz2pow.InputDelay = 1
frfBode(G_uz2pow, freqs, F1, 'Hz', '--m');

%%
R2 = 1.7e6;
R1 = 28.4e6;
Vdiv = R2/(R1 + R2);
Vdiv_gain = 1/Vdiv % from resistor measurement
Vdiv_gain = 8.196/0.451; % measured from labview + 9v battery
G_uz2pow_rescaled = G_uz2pow*Vdiv_gain;

modelFit.models.G_uz2pow = G_uz2pow;
modelFit.models.G_uz2pow_rescaled = G_uz2pow_rescaled;
modelFit.models.Vdiv_gain = Vdiv_gain;

save(modelFit_file, 'modelFit')
%% 

C1 = 4e-6;
Imax = 100e-3;
dVmax_high = (Ts/C1)*Imax
dVmax_low = dVmax_high/dcgain(G_uz2pow_rescaled)
%%
% Now, validate in time time domain. The point here is to make sure that
% our model and scaling makes sense, not to explore the limits. Do that in
% measure_slew_limit.m
ms = 1e3;

To_half = 40*Ts;
To = 2*To_half;
wo = 2*pi/To;

t0 = [0:1:(3*To/Ts)-1]'*Ts;
du = 0.15;
du_vec = sign(sin(wo*t0))*du;
u_vec = [cumsum(du_vec); zeros(800, 1)];
t = (0:length(u_vec)-1)'*Ts;



F4 = figure(4); clf
subplot(4,1,1)
plot(t*ms, u_vec)
grid on
title('(LV) u(k)', 'interpreter', 'latex')
xlabel('time [ms]')

subplot(4,1,2)
plot(t(1:end-1)*ms, diff(u_vec))
grid on
title('(LV) $\Delta u(k)$', 'interpreter', 'latex')
xlabel('time [ms]')

% du_pow = lsim(Gpow/dcgain(Gpow), diff(u), t(1:end-1));
du_pow = lsim(G_uz2pow_rescaled, diff(u_vec), t(1:end-1));
subplot(4,1,3)
plot(t(1:end-1)*ms, du_pow)
grid on
xlabel('time [ms]')

subplot(4,1,4)
plot(t(1:end-1)*ms, cumsum(du_pow))
grid on;
xlabel('time [ms]')
ylabel('$y_{pow}$')
hold on
%%


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
    vipath = 'C:\Users\arnold\Documents\MATLAB\afm_mpc_journal\labview\play_nPoint_id_slew_OL_FIFO.vi'

    [e, vi] = setupVI(vipath, 'Abort', 0,...
                'umax', 6, 'data_out_path', slewfpath_out,...
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
upow_exp = upow_exp - upow_exp(1); 

% upow_exp = data(:,3)*Vdiv_gain/dcgain(Gpow);
texp = [0:1:length(u_exp)-1]'*Ts;

figure(F4);
subplot(4,1,1), hold on;
plot(texp*ms, u_exp, '--k')


subplot(4,1,3), hold on
plot(texp(1:end-1)*ms, diff(upow_exp), '--')
title('(HV) $\Delta u(k)$', 'interpreter', 'latex')


subplot(4,1,4)
plot(texp*ms, upow_exp, '--')
title('(HV) $y_{pow}$', 'interpreter', 'latex')
grid on





