clear, clc

saveon = false;

addpath(fullfile(getMatPath, 'afm_mpc_journal', 'functions'));

data_root = fullfile(PATHS.exp, 'step-exps/many_steps_data_rand_27-Sep-2018_01');

% ls(data_root)
% 
% names_sim_v1 = {'single_step_linfxp_sim_const-sig-same-sig_08-30-2018.mat',...
% 'single_step_mpcfxp_sim_const-sig-same-sig_08-30-2018.mat'};
TOL = 0.01;
TOL = 14/512;


names_exp_v1 = {'single_step_lin_EXP_const-sig-same-sig_09-27-201801.mat',...
'single_step_mpc_EXP_const-sig-same-sig_09-27-201801.mat'};

load(fullfile(data_root, names_exp_v1{1}));

lin_exp_v1 = afm_exp_lin;
% lin_exp_v1.controller;

load(fullfile(data_root, names_exp_v1{2}));
mpc_exp_v1 = afm_exp_mpc;

mpc_exp_v1.Color = 'r';
mpc_exp_v1.LineStyle = '--';
mpc_exp_v1.yscaling = 5;
mpc_exp_v1.yunits = '$\mu$m';

tsmpc1 = mpc_exp_v1.settle_time(TOL, 'abs', 0)*1000;
tslin1  = lin_exp_v1.settle_time(TOL, 'abs', 0)*1000;


mpc_exp_v1.name = sprintf('MPC-CS, $\\gamma=0.11$, settle-time=%.2f [ms]', tsmpc1);
lin_exp_v1.name = sprintf('SLF-CS, $\\gamma=0.11$, settle-time=%.2f [ms]', tslin1);
lin_exp_v1.yscaling = 5;
lin_exp_v1.yunits = '$\mu$m';


names_sim_v2 = {'single_step_linfxp_sim_const-sig-same-sig_09-27-201805.mat',...
'single_step_mpcfxp_sim_const-sig-same-sig_09-27-201805.mat'};

names_exp_v2 = {'single_step_lin_EXP_const-sig-same-sig_09-27-201802.mat',...
  'single_step_mpc_EXP_const-sig-same-sig_09-27-201802.mat'};
load(fullfile(data_root, names_exp_v2{1}))
load(fullfile(data_root, names_exp_v2{2}))


lin_exp_v2 = afm_exp_lin;
mpc_exp_v2 = afm_exp_mpc;

figure(1); hold on;
ax1 = gca();

hlin2 = lin_exp_v2.ploty(ax1);
hlin2 = mpc_exp_v2.ploty(ax1);
mpc_exp_v2.step_ref.plot_settle_boundary(ax1, TOL, 'abs')

%%

tslin2  = lin_exp_v2.settle_time(TOL, 'abs', 0)*1000;

lin_exp_v2.name = sprintf('SLF-CS, $\\gamma=12$, settle-time=%.2f [ms]', tslin2);
lin_exp_v2.Color = 'k';
lin_exp_v2.LineStyle = '-';
lin_exp_v2.yscaling = 5;
lin_exp_v2.yunits = '$\mu$m';

%
width = 3.45;
height = 3;
F1 = mkfig(1, width, height); clf;
lft = 0.1027;
bt1 = 0.5638;
bt2 = 0.1100;
ht1 = 0.3953;
wd = 0.8610;
ax1 = axes('Position', [lft, bt1, wd, ht1+.02]);
ax2 = axes('Position', [lft, bt2, wd, ht1]);
% subplot(2,1,1);
% ax2 = subplot(2,1,2);

mpc_exp_v1.Color = 'r';
mpc_exp_v1.LineStyle = '-';
lin_exp_v1.Color = 'b';
lin_exp_v1.LineStyle = '-';
lin_exp_v2.Color = 'k';
lin_exp_v2.LineStyle = '--';


hmpc1 = mpc_exp_v1.ploty(ax1);
hlin1 = lin_exp_v1.ploty(ax1);
hlin2 = lin_exp_v2.ploty(ax1);

grid(ax1, 'on')

mpc_exp_v1.step_ref.yscaling = 5;
mpc_exp_v1.step_ref.plot_settle_boundary(ax1, TOL, 'abs')

lin_exp_v1.plotdu(ax2);
mpc_exp_v1.plotdu(ax2);
lin_exp_v2.plotdu(ax2);

grid(ax2, 'on')

tf = 0.015;
tstart = 0.00;
xlim(ax1, [tstart, tf])
xlim(ax2, [tstart, tf])

ylim(ax1, [1, 1.5]*5)

leg1 = legend([hmpc1, hlin1, hlin2]);
set(leg1, 'Position', [0.2664 0.5736 0.7122 0.1449], 'Box', 'off')
ylabel(ax1, '$y$ [$\mu$m]')
ylabel(ax2, '$\Delta u_k$')
xlabel(ax2, 'time [s]')

set(ax1, 'XTickLabel', [])
if saveon
  saveas(F1, fullfile(PATHS.jfig, 'disingenuous.svg'))
end
%%
% now simulation

names_sim_v1 = {'single_step_linfxp_sim_const-sig-same-sig_08-30-2018.mat',...
'single_step_mpcfxp_sim_const-sig-same-sig_08-30-2018.mat'};


load(fullfile(data_root, names_sim_v1{1}));
lin_sim_v1 = sim_exp_fxpl;
lin_sim_v1.controller;

load(fullfile(data_root, names_sim_v1{2}));

mpc_sim_v1 = sim_exp_fxpm;

mpc_sim_v1.Color = 'r';
mpc_sim_v1.LineStyle = '--';
mpc_sim_v1.yscaling = 5;
mpc_sim_v1.yunits = '$\mu$m';

tsmpc1 = mpc_sim_v1.settle_time(TOL, 'abs', 0)*1000
tslin1  = lin_sim_v1.settle_time(TOL, 'abs', 0)*1000


mpc_sim_v1.name = sprintf('MPC-CS, $\\gamma=0.11$, settle-time=%.2f', tsmpc1)
lin_sim_v1.name = sprintf('SLF-CS, $\\gamma=0.11$, settle-time=%.2f', tslin1)
lin_sim_v1.yscaling = 5;
lin_sim_v1.yunits = '$\mu$m';

lin_sim_v1.Color = 'b';
lin_sim_v1.LineStyle = '-';

names_sim_v2 = {'single_step_linfxp_sim_const-sig-same-sig_08-30-201802.mat',...
'single_step_mpcfxp_sim_const-sig-same-sig_08-30-201802.mat'};

names_sim_v2 = {'single_step_linfxp_sim_const-sig-same-sig_08-30-201802.mat'};
load(fullfile(data_root, names_sim_v2{1}))
lin_sim_v2 = sim_exp_fxpl;
tslin2  = lin_sim_v2.settle_time(TOL, 'abs', 0)*1000

lin_sim_v2.name = sprintf('SLF-CS, $\\gamma=12$, settle-time=%.2f', tslin2)
lin_sim_v2.Color = 'k';
lin_sim_v2.LineStyle = '-';
lin_sim_v2.yscaling = 5;
lin_sim_v2.yunits = '$\mu$m';


width = 3.45;
height = 3;
F2 = mkfig(2, width, height); clf;
lft = 0.1027;
bt1 = 0.5638;
bt2 = 0.1100;
ht1 = 0.3953;
wd = 0.8610;
ax1 = axes('Position', [lft, bt1, wd, ht1+.02]);
ax2 = axes('Position', [lft, bt2, wd, ht1]);

hlin2 = lin_sim_v2.ploty(ax1);
hmpc1 = mpc_sim_v1.ploty(ax1);
hlin1 = lin_sim_v1.ploty(ax1);

grid(ax1, 'on')

mpc_sim_v1.step_ref.yscaling = 5;
mpc_sim_v1.step_ref.plot_settle_boundary(ax1, TOL, 'abs')

lin_sim_v1.plotdu(ax2);
lin_sim_v2.plotdu(ax2);
mpc_sim_v1.plotdu(ax2);

grid(ax2, 'on')

tf = 0.015;
xlim(ax1, [0, tf])
xlim(ax2, [0, tf])

ylim(ax1, [1, 1.5]*5)

leg2 = legend([hmpc1, hlin1, hlin2]);
set(leg2, 'Position', [0.2513 0.6500 0.7122 0.1449], 'Box', 'off')
ylabel(ax1, '$y$ [$\mu$m]')
ylabel(ax2, '$\Delta u_k$')
xlabel(ax2, 'time [s]')

set(ax1, 'XTickLabel', [])


%%
saveas(F1, fullfile(PATHS.jfig, 'disingenuous_exp.svg'))
saveas(F2, fullfile(PATHS.jfig, 'disingenuous_sim.svg'))
