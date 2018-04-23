clear
clc

modelFit_file = 'C:\Users\arnold\Documents\labview\sysID\data\x-axis_sines_info_intsamps_amp_1p0_out_4-2-2018offs_0-01.mat';
load(modelFit_file)

Ts = modelFit.frf.Ts;
G_uz2pow_rescaled = modelFit.models.G_uz2pow_rescaled;


ms = 1e3;

%% Case 1:
% The signal is not too aggressive. We expect nearly perfect match between
% power amp model and output. 
clc
To_half = 28*Ts;
To = 2*To_half;
wo = 2*pi/To;

Ts = 40e-6;

n_space = 15000;
% Define a sequence of impulses. Their integral is a sequence of steps.
N_imp = 8;
if 1
N_imp_vec= (n_space:n_space:N_imp*n_space)';
u_vec = zeros(N_imp*n_space, 1);
u_vec(N_imp_vec(1:N_imp/2)) = 1;
u_vec(N_imp_vec(N_imp/2+1:end)) = -1;
u_vec = repmat(cumsum(u_vec), 3,1);

t_vec = (0:length(u_vec)-1)'*Ts;
plot(t_vec, u_vec);
else
    [rast] = raster(1/(2*n_space*Ts), Ts, N_imp*n_space*Ts)
    u_vec = rast.Data*2;
    t_vec = rast.Time;
end

%
clc

F4 = figure(4); clf
subplot(3,1,1)
plot(t_vec*ms, u_vec)
grid on
title('(LV) u(k)', 'interpreter', 'latex')
xlabel('time [ms]')

subplot(3,1,2)
plot(t_vec(1:end-1)*ms, diff(u_vec))
grid on
title('(LV) $\Delta u(k)$', 'interpreter', 'latex')
xlabel('time [ms]')



slewfname_in = 'hyst_id_datain_4-20-2018_02.csv';
% slewfpath_in = fullfile(pwd, slewfname_in)
slewfpath_in = fullfile(PATHS.exp, 'sysID', slewfname_in);

slewfname_out = 'hyst_id_4-20-2018_02.csv';
slewfpath_out = fullfile(PATHS.exp, 'sysID', slewfname_out);
% slewfpath_out = fullfile(pwd, slewfname_out);
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
clc
figure(2); clf

dat = csvread(slewfpath_out);

u_exp = dat(:,1);
yx_exp = dat(:,2) - dat(1,2);
Vpow_exp = dat(:,3);
I_exp = dat(:,4);
t_exp = (0:length(u_exp)-1)'*Ts;

plot(t_exp, u_exp*.68)
hold on
plot(t_exp, yx_exp)
grid on

figure(3), clf
plot(u_exp, yx_exp, 'o')
hold on
plot(u_exp(1), yx_exp(1), 'o')
grid on
xlabel('u')
ylabel('y')


