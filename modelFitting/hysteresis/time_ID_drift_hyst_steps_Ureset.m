clear
clc

modelFit_file = fullfile(PATHS.sysid, 'FRF_data_current_stage2.mat');
load(modelFit_file)

Ts = modelFit.frf.Ts;
G_uz2stage = modelFit.models.G_uz2stage;


ms = 1e3;

%Case 1:

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
%%
n_space = 15000;
% Define a sequence of impulses. Their integral is a sequence of steps.
figure(1); clf; 
N_imp = 8;
if 0
N_imp_vec= (n_space:n_space:N_imp*n_space)';
u_vec = zeros(N_imp*n_space, 1);
u_vec(N_imp_vec(1:N_imp/2)) = 1;
u_vec(N_imp_vec(N_imp/2+1:end)) = -1;
u_vec = repmat(cumsum(u_vec), 3,1);

t_vec = (0:length(u_vec)-1)'*Ts;
plot(t_vec, u_vec);
grid on
else
    [rast] = raster(1/(2*n_space*Ts), Ts, N_imp*n_space*Ts);
    u_vec = rast.Data*umax;
    t_vec = rast.Time;
    plot(t_vec, u_vec);
    grid on;
end

u_vec_both = [u_reset; u_vec];
t_vec_both = (0:length(u_vec_both)-1)'*Ts;
figure(2); clf
plot(t_vec_both, u_vec_both)
grid on

figure(3)
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
    csvwrite(slewfpath_in, u_vec_both);
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
yx_exp = dat(:,2); % - dat(1,2);
Vpow_exp = dat(:,3);
I_exp = dat(:,4);
t_exp = (0:length(u_exp)-1)'*Ts;

plot(t_exp, u_exp*.68)
hold on
plot(t_exp, yx_exp)
grid on

idx_ureset_end = length(u_reset);
figure(3), clf
yx = yx_exp(idx_ureset_end+1:end) - yx_exp(idx_ureset_end+1);
plot(u_exp(idx_ureset_end+1:end), yx, 'o')
hold on
plot(u_exp(1), yx_exp(1), 'o')
grid on
xlabel('u')
ylabel('y')

%%
clc

n = 7;
r = linspace(0, umax, n);
w = ones(n, 1);


hyst_fun = @(w, y0) PIHyst.hyst_play_op(u_exp(idx_ureset_end+1:end), r, w, y0);
cost_fun = @(w) hyst_fun(w, w*0) - yx;
w_est = lsqnonlin(cost_fun, w);

y_est = PIHyst.hyst_play_op(u_exp(idx_ureset_end+1:end), r, w_est, w*0);

figure(5)
plot(yx)
hold on, grid on
plot(y_est, '--')

figure(3)
plot(u_exp(idx_ureset_end+1:end), y_est, '--')

%%

N_imp = 8;

N_imp_vec= (n_space:n_space:N_imp*n_space)';
u_vec = zeros(N_imp*n_space, 1);
u_vec(N_imp_vec(1:N_imp/2)) = 1;
u_vec(N_imp_vec(N_imp/2+1:end)) = -1;
u_vec = repmat(cumsum(u_vec), 3,1);

t_vec = (0:length(u_vec)-1)'*Ts;

figure(6)
plot(t_vec, u_vec);
grid on

u_vec_both = [u_reset; u_vec];
t_vec_both = (0:length(u_vec_both)-1)'*Ts;
figure(2); clf
plot(t_vec_both, u_vec_both)
grid on

slewfname_in = 'hyst_id_datain_4-20-2018_03.csv';
slewfpath_in = fullfile(PATHS.exp, 'sysID', slewfname_in);

slewfname_out = 'hyst_id_4-20-2018_03.csv';
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
%%
dat_steps = csvread(slewfpath_out);

u_exp_steps = dat_steps(:,1);
yx_exp_steps = dat_steps(:,2); 
t_exp = (0:length(u_exp)-1)'*Ts;

idx_ureset_end = length(u_reset);

ux_steps = u_exp_steps(idx_ureset_end+1:end);
yx_steps = yx_exp_steps(idx_ureset_end+1:end) - yx_exp_steps(idx_ureset_end+1);

figure(4), clf
plot(yx_steps)
hold on
grid on
xlabel('k')
ylabel('y')

y_est_steps = PIHyst.hyst_play_op(ux_steps, r, w_est, w*0);

hyst_fun = @(w, y0) PIHyst.hyst_play_op(ux_steps, r, w, y0);
cost_fun = @(w) hyst_fun(w, w*0) - yx_steps;
w_est_steps = lsqnonlin(cost_fun, w_est);
y_est_steps_custom = PIHyst.hyst_play_op(ux_steps, r, w_est_steps, w*0);

plot(y_est_steps, '--')
plot(y_est_steps_custom, '-g')

% suppose we have a third order drift model. Lets first fit a model to the
% first step, and use that as the initial guess to optimize both the hysteresis
% model and drift model together.

%%
idx_steps = find(diff(ux_steps) ~=0);
idx_step1_end = idx_steps(2)-1;

u1 = ux_steps(1:idx_step1_end);
y1 = yx_steps(1:idx_step1_end);
t1_vec = (0:length(y1)-1)'*Ts;

figure(20); clf
plot(y1)
hold on
plot(u1)

% Imagine that the output of our vibrational model is the input to the drift
% model. This lets us separate the affects.
% [Gvib, gdrift1] = eject_gdrift(G_uz2stage);

% u_interm = lsim(Gvib, u1, t1_vec);
% 
% plot(u_interm)

theta0 = [0.99928,0.99828 0.99840,0.99940 1.2];
% from previous run:
theta0 = [0.9922    0.9997    0.9997    0.9927    1.0152];
np = 2;

gdrift_cost = @(theta)fit_gdrift(theta, G_uz2stage, y1, u1, t1_vec, np);

opts = optimoptions(@lsqnonlin);
opts.MaxFunctionEvaluations = 5000;
theta = lsqnonlin(gdrift_cost, theta0, [], [], opts);
%
gdrift2 = zpk(theta(np+1:end-1), theta(1:np), theta(end), Ts);

ydrift_est0 = lsim(G_uz2stage*gdrift2, u1, t1_vec);

plot(ydrift_est0, '--r')

%%
% Now, we want to use the initial system, and optimize (simultaneously) the
% hysteresis and the drift.
clc
t_vec_steps = (0:length(yx_steps)-1)'*Ts;

theta0_both = [theta'; w_est];
idx_theta_drft = (1:2*np+1)';
idx_theta_hyst = (2*np+2:length(theta0_both))';

%%
G = G_uz2stage*gdrift2;
% G = G/dcgain(G);

u_hyst_fun = @(theta) PIHyst.hyst_play_op(ux_steps, r, theta, theta*0);

cost_hyst_only = @(theta) lsim(G, u_hyst_fun(theta),... 
                t_vec_steps) - yx_steps; 
opts = optimoptions(@lsqnonlin);
opts.MaxFunctionEvaluations = 1000;
opts.Display = 'iter';
theta_hyst3 = lsqnonlin(cost_hyst_only, w_est_steps, [],[], opts);
%%
u_est_steps_all = PIHyst.hyst_play_op(ux_steps, r, theta_hyst3, w*0);

y_est_all = lsim(G, u_est_steps_all, t_vec_steps);
%
figure(15); clf;
plot(t_vec_steps, yx_steps)
hold on
plot(t_vec_steps, y_est_all, 'r')
plot(t_vec_steps, ux_steps*dcgain(G_uz2stage))
grid on

%%
gdrift_cost_all = @(theta, u)fit_gdrift(theta, G_uz2stage, yx_steps, u, t_vec_steps, np);

hyst_fun = @(w) PIHyst.hyst_play_op(ux_steps, r, w, w*0);

cost_fun = @(theta) gdrift_cost_all(theta(idx_theta_drft), hyst_fun(theta(idx_theta_hyst)));

opts = optimoptions(@lsqnonlin);
opts.MaxFunctionEvaluations = 1000;
opts.Display = 'iter';
theta_est_steps = lsqnonlin(cost_fun, theta0_both,[],[], opts);
%%
theta_est_steps = lsqnonlin(cost_fun, theta_est_steps,[],[], opts);
%%
theta_drft = theta_est_steps(idx_theta_drft);
theta_hyst = theta_est_steps(idx_theta_hyst);
gdrift_all = zpk(theta_drft(np+1:end-1), theta_drft(1:np), theta_drft(end), Ts)


cost_fun2 = @(theta) gdrift_cost_all(theta(idx_theta_drft), hyst_fun(theta_hyst));

theta_drft2 = lsqnonlin(cost_fun2, theta_drft,[],[], opts);

%%
gdrift_all = zpk(theta_drft2(np+1:end-1), theta_drft2(1:np), theta_drft2(end), Ts)
%%
% fun = @(theta) gdrift_cost_all(theta(idx_theta_drft), hyst_fun(theta(idx_theta_hyst)))

% y_est_all = fun(theta_est_steps);

u_est_steps_all = PIHyst.hyst_play_op(ux_steps, r, theta_hyst, w*0);

y_est_all = lsim(G_uz2stage*gdrift_all, u_est_steps_all, t_vec_steps);
%
figure(15); clf;
plot(t_vec_steps, yx_steps)
hold on
plot(t_vec_steps, y_est_all)
plot(t_vec_steps, ux_steps*dcgain(G_uz2stage))
grid on


% % %%
% % [rp, wp] = PIHyst.invert_hyst_PI(r, w_est_steps);
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



