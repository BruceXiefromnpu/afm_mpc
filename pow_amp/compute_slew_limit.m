clear
clc


Imax = 100e-3; % 100mA
Ts = 40e-6;
% assert(1/Ts == 25e3);
C = 3.8e-6; % muF

Kamp = 9;  % 20v range to 180 v range

del_Vhigh_max = (Ts/C)*Imax

del_Vlow_max = del_Vhigh_max/Kamp


%%

Vdiv = 2.56/(19.61+2.56);
Vdiv_gain = 1/Vdiv;
load( fullfile(PATHS.exp, 'x-axis_sines_info_out_2-8-2018-01.mat'))
Gpow = modelFit.models.G_uz2pow*Vdiv_gain;
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


du_pow = lsim(Gpow, diff(u), t(1:end-1));
subplot(4,1,3)
plot(t(1:end-1)*ms, du_pow)
grid on
%%
slewfname_in = sprintf('data/slewexp_datain_%0.2f.csv', max(diff(triang.Data)));
slewfpath_in = fullfile(PATHS.MPCJ_root, slewfname_in);

slewfname_out = sprintf('data/slewexp_dataout_%0.2f.csv', max(diff(triang.Data)));
slewfpath_out = fullfile(PATHS.MPCJ_root, slewfname_out);



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
%%
data = csvread(slewfpath_out);
size(data)

u_exp = data(:,1);
y_exp = data(:,2);
upow_exp = data(:,3)*Vdiv_gain;
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

figure(5);
plot(y_exp)






