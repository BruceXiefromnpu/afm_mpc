clear
clc

modelFit_file = fullfile(PATHS.sysid, 'FRF_data_current_stage2.mat');
load(modelFit_file)

Ts = modelFit.frf.Ts;
G_uz2stage = modelFit.models.G_uz2stage;


ms = 1e3;

% Case 1:

clc
To_half = 28*Ts;
To = 2*To_half;
wo = 2*pi/To;

% Build the u-reset.
t1 = 18;
tf = 100;
umax = 5;
k1 = 0.15;
u_reset = PIHyst.gen_reset_u(t1, tf, Ts, k1, umax);
figure(100)
plot(u_reset)
%
n_step_start = floor(20/Ts);
n_step = floor(20/Ts);
step_amp = 0.25;
u_vec = [zeros(n_step_start, 1);
        step_amp*ones(n_step, 1)];
t_vec = (0:length(u_vec)-1)'*Ts;

u_vec_both = [u_reset; u_vec];
t_vec_both = (0:length(u_vec_both)-1)'*Ts;

figure(2); clf
plot(t_vec_both, u_vec_both)
grid on


slewfname_in = 'drift_id_datain_4-23-2018_02.csv';
slewfpath_in = fullfile(PATHS.exp, 'sysID', slewfname_in);

slewfname_out = 'drift_id_4-23-2018_02.csv';
slewfpath_out = fullfile(PATHS.exp, 'sysID', slewfname_out);
%%
% -----------------------RUN THE Experiment--------------------------------
if 1
    csvwrite(slewfpath_in, u_vec_both);
    clear vi;
    vipath = 'C:\Users\arnold\Documents\MATLAB\afm_mpc_journal\labview\play_nPoint_id_slew_OL_FIFO.vi'

    [e, vi] = setupVI(vipath, 'Abort', 0,...
                'umax', 9, 'data_out_path', slewfpath_out,...
                'traj_in_path', slewfpath_in, 'TsTicks', 1600);
    vi.Run
end
%
clc
figure(2); clf

dat = csvread(slewfpath_out);

u_exp = dat(:,1);
yx_exp = dat(:,2); % - dat(1,2);
Vpow_exp = dat(:,3);
I_exp = dat(:,4);
t_exp = (0:length(u_exp)-1)'*Ts;

plot(t_exp, u_exp*dcgain(G_uz2stage))
hold on
plot(t_exp, yx_exp)
grid on

idx_ureset_end = length(u_reset);
figure(3), clf
yx = yx_exp(idx_ureset_end+1:end) - yx_exp(idx_ureset_end+1);


%%



