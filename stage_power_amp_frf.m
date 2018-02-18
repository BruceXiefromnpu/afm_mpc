clc
load('C:\Users\arnold\Documents\labview\sysID\data\x-axis_sines_info_out_12-10-2017-04.mat');

mf_old = modelFit;

load('C:\Users\arnold\Documents\labview\sysID\data\x-axis_sines_info_out_2-8-2018-01.mat');
mf_new = modelFit;



F1 = figure(1);

plt_hands1 = frfBode(mf_new.frf.Gx_ux_OL_frf, F1, '-r', 0);
plt_hands2 = frfBode(mf_old.frf.Gx_ux_OL_frf, F1, '--k', 0);
plt_hands3 = frfBode(squeeze(mf_new.frf.G_pow_frf), mf_new.frf.freq_s, F1, '-b', 'Hz');

hand_mag_stage_new = plt_hands1.h_mag;
hand_mag_stage_new.DisplayName = 'stage + amp: new';

hand_mag_stage_old = plt_hands2.h_mag;
hand_mag_stage_old.DisplayName = 'stage + amp: old';

hand_mag_amp_new = plt_hands3;
hand_mag_amp_new.DisplayName = 'power amp';

leg1 = legend([hand_mag_stage_new, hand_mag_stage_old, hand_mag_amp_new]);
set(leg1, 'Location', 'SouthWest')

ax1 = plt_hands1.AX(1);
F1.CurrentAxes = ax1;

set(ax1, 'YTick', [-100:20:0])
%%

F2 = figure(2);
[frf_pow, freqs]= monotonicFRF(squeeze(mf_new.frf.G_pow_frf), mf_new.frf.freq_s);
frf_stage = squeeze(mf_new.frf.Gx_ux_OL_frf.ResponseData)./frf_pow;

frfBode(frf_stage, freqs, F2, 'r', 'Hz')


%%
pow_frd = frd(frf_pow, freqs*2*pi, Ts);

Gpow_est = tfest(pow_frd, 2, 0, 'InputDelay', 5, 'Ts', Ts);
%%
clc
Gpow_est_frf = squeeze(freqresp(Gpow_est, freqs*2*pi));


F3 = figure(3);
plt_hands3 = frfBode(squeeze(mf_new.frf.G_pow_frf), mf_new.frf.freq_s, F3, '-b', 'Hz');
frfBode(Gpow_est_frf, freqs, F3, '--k', 'Hz');

save('gpow_tf.mat', 'Gpow_est')

