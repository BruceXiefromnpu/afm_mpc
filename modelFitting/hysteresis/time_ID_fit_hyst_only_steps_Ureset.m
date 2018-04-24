% This script tries to fit ONly the hysteresis model. It uses the separeted
% Gvib and drift model fit in fit_drift.m

clear
clc

load('drift_data.mat')

Ts = G_uz2stage.Ts;

ms = 1e3;


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





% Load the steps data with the initializing reset.
slewfname_in = 'hyst_id_datain_4-20-2018_03.csv';
slewfpath_in = fullfile(PATHS.exp, 'sysID', slewfname_in);

slewfname_out = 'hyst_id_4-20-2018_03.csv';
slewfpath_out = fullfile(PATHS.exp, 'sysID', slewfname_out);


dat = csvread(slewfpath_out);
%%
u_exp = dat(:,1);
u_exp(end-3:end) = u_exp(end-4);
y_exp = dat(:,2); % - dat(1,2);
Vpow_exp = dat(:,3);
I_exp = dat(:,4);
t_exp = (0:length(u_exp)-1)'*Ts;

figure(2); clf

plot(t_exp, u_exp*.68)
hold on
plot(t_exp, y_exp)
grid on

k_start = length(u_reset);

ux = u_exp(k_start:end);
yx = y_exp(k_start:end);
k1 = find(diff(ux) ~=0, 1, 'first') - 1;
yx = yx - mean(yx(1:k1));
tvec = (0:length(yx)-1)'*Ts;

figure(3); clf
plot(tvec, yx)
hold on
plot(tvec, ux)


%%


n = 7;
r = linspace(0, umax, n);
w = ones(n, 1);

% hyst_fun = @(w, y0) PIHyst.hyst_play_op(u_exp(idx_ureset_end+1:end), r, w, y0);

G = Gvib*gdrift;

% u_hyst_fun = @(theta) PIHyst.hyst_play_op(ux, r, theta(1:n), theta(n+1:end));
u_hyst_fun = @(theta) PIHyst.hyst_play_op(ux, r, theta(1:n), w*0);

cost_hyst_only = @(theta) lsim(G, u_hyst_fun(theta),... 
                tvec) - yx; 
opts = optimoptions(@lsqnonlin);
opts.MaxFunctionEvaluations = 1000;
opts.Display = 'iter';
theta_hyst = lsqnonlin(cost_hyst_only, w, [],[], opts);

% u_est_steps_all = PIHyst.hyst_play_op(ux, r, theta_hyst(1:n), theta_hyst(n+1:end));
u_est_steps_all = PIHyst.hyst_play_op(ux, r, theta_hyst(1:n), w*0);

y_est_all = lsim(G, u_est_steps_all, tvec);
%%
figure(15); clf;
h1 = plot(tvec, yx);
hold on
h2 = plot(tvec, y_est_all, 'r');
h3 = plot(tvec, ux*dcgain(G), 'LineWidth', 1.5);

h1.DisplayName = 'Experimental';
h2.DisplayName = 'Fit';
h3.DisplayName = 'u*dcgain(G)';
legend([h1, h2, h3])
grid on




%%
[rp, wp] = PIHyst.invert_hyst_PI(r, theta_hyst);

save('steps_hyst_model.mat', 'rp', 'wp', 'r', 'theta_hyst', 'umax', 'Gvib', 'gdrift')



% % 
% % u_inv = PIHyst.hyst_play_op(u_vec, rp, wp, wp*0);
% % u_vec_both = [u_reset; u_inv];
% % 
% % slewfname_in = 'hyst_id_datain_4-20-2018_04.csv';
% % slewfpath_in = fullfile(PATHS.exp, 'sysID', slewfname_in);
% % 
% % slewfname_out = 'hyst_id_4-20-2018_04.csv';
% % slewfpath_out = fullfile(PATHS.exp, 'sysID', slewfname_out);
% % 
% % %%
% % % -----------------------RUN THE Experiment--------------------------------
% % if 1
% %     csvwrite(slewfpath_in, u_vec_both);
% %     clear vi;
% %     vipath = 'C:\Users\arnold\Documents\MATLAB\afm_mpc_journal\labview\play_nPoint_id_slew_OL_FIFO.vi'
% % 
% %     [e, vi] = setupVI(vipath, 'Abort', 0,...
% %                 'umax', 9, 'data_out_path', slewfpath_out,...
% %                 'traj_in_path', slewfpath_in, 'TsTicks', 1600);
% %     vi.Run
% % end
% % %%
% % dat_steps = csvread(slewfpath_out);
% % 
% % u_exp_steps = dat_steps(:,1);
% % yx_exp_steps = dat_steps(:,2); 
% % t_exp = (0:length(u_exp)-1)'*Ts;
% % 
% % idx_ureset_end = length(u_reset);
% % %%
% % ux_steps = u_exp_steps(idx_ureset_end+1:end);
% % yx_steps = yx_exp_steps(idx_ureset_end+1:end) - yx_exp_steps(idx_ureset_end+1);
% % 
% % figure(5), clf
% % plot(yx_steps)
% % hold on
% % grid on
% % plot([u_vec, ux_steps(1:end-1)], 'LineWidth', 1.5)
% % xlabel('k')
% % ylabel('y')



