

clear, clc


C_stage = 3.8e-6;
Imax = 100e-3;
Ts = 40e-6;


Vdiv = 2.56/(19.61+2.56);
Vdiv_gain = 1/Vdiv; % from resistor measurement
Vdiv_gain = 10.6; % From 9v battery measurement
load( fullfile(PATHS.exp, 'x-axis_sines_info_out_2-8-2018-01.mat'))
Gpow_hat = modelFit.models.G_uz2pow
Gpow = Gpow_hat*Vdiv_gain;
% dcgain(Gpow)

dVmax = (Ts/C_stage)*Imax
% dVmax = (Ts/(C_stage*dcgain(Gpow)))*Imax


%%
Ts = 40e-6;
ms = 1e3;
period = (60*Ts)
triang = raster(1/period, Ts, period);
triang.Data = triang.Data*6;
u = [triang.Data; zeros(400,1)];
t = [0:1:length(u)-1]'*Ts;

F4 = figure(4); clf
subplot(4,1,1)
plot(t*ms, u)
grid on
title('(LV) u(k)', 'interpreter', 'latex')

subplot(4,1,2)
plot(t(1:end-1)*ms, diff(u))
grid on
title('(LV) $\Delta u(k)$', 'interpreter', 'latex')

% du_pow = lsim(Gpow/dcgain(Gpow), diff(u), t(1:end-1));
du_pow = lsim(Gpow, diff(u), t(1:end-1));
subplot(4,1,3)
plot(t(1:end-1)*ms, du_pow)
grid on
%%
slewfname_in = sprintf('slewexp_datain_%0.2f.csv', max(diff(triang.Data)));
slewfpath_in = fullfile(PATHS.exp, 'sysID', slewfname_in);

slewfname_out = sprintf('slewexp_dataout_%0.2f.csv', max(diff(triang.Data)));
slewfpath_out = fullfile(PATHS.exp, 'sysID', slewfname_out);

%
% -----------------------RUN THE Experiment--------------------------------
if 0
    csvwrite(slewfname_in, u);
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
% upow_exp = data(:,3)*Vdiv_gain/dcgain(Gpow);
texp = [0:1:length(u_exp)-1]'*Ts;

figure(F4);
subplot(4,1,1), hold on;
plot(texp*ms, u_exp, '--')


subplot(4,1,3), hold on
plot(texp(1:end-1)*ms, diff(upow_exp), '--')
title('(HV) $\Delta u(k)$', 'interpreter', 'latex')


subplot(4,1,4)
plot(texp*ms, upow_exp)
title('(HV) $ u_{lpf}(k)$', 'interpreter', 'latex')
grid on

%
% Now plot a nice figure 
clc
F6 = mkfig(6, 5, 4);
clf

% subplot(1,2,1)
% ax1 = gca;
lft1 = 0.08, 
lft2 = 0.5;
bt = 0.11;
wd1 = 0.3347;
wd2 = 0.45;
ht = 0.815;

ax1 = axes('Position', [lft1, bt, wd1, ht])
plot(t*ms, u, 'LineWidth', 2);
grid on
title('(LV) u(k)', 'interpreter', 'latex')
xlim([0, 5])
xlabel('time [ms]')


% subplot(1,2,2)
% ax2 = gca;
ax2 = axes('Position', [lft2, bt, wd2, ht])

hold on, grid on;
h1 = plot(t(1:end-1)*ms, du_pow, '-');

h2 = plot(texp(1:end-1)*ms, diff(upow_exp), '--', 'LineWidth', 2);
xlm = xlim;
plot(xlm, [dVmax, dVmax], ':k', 'LineWidth', 2);
h3 = plot(xlm, -[dVmax, dVmax], ':k', 'LineWidth', 2);

title('(HV) $\Delta u(k)$', 'interpreter', 'latex')
xlim([0, 5])
xlabel('time [ms]')

h1.DisplayName = '$\Delta u(k)$ (model)';
h2.DisplayName = '$\Delta u(k)$ (measured)';
h3.DisplayName = 'Predicted $\Delta u_{max}$';

leg1 = legend(gca, [h1, h2, h3]);
leg1.Position = [0.7043 0.7642 0.2846 0.1255];
saveon = 1;
if saveon
    saveas(F6, fullfile(PATHS.jfig, 'dumax_measure.svg'));
end




