% Plot the stage vs simple cap circuit.
clear
clc

MF_stage = load('FRF_data_current_stage2.mat');
MF_cap   = load('FRF_data_current_caponly2.mat');

MF_stage = MF_stage.modelFit;
MF_cap = MF_cap.modelFit;

freqs_cap = MF_cap.frf.freq_s;
freqs_stage = MF_stage.frf.freq_s;

F1 = figure(10);
F2 = figure(20);


h1 = frfBode(MF_stage.frf.G_powV2powI, freqs_stage, F1, '-r', 'Hz');
h1.DisplayName = 'powV2powI, stage';

h2 = frfBode(MF_cap.frf.G_powV2powI, freqs_cap, F1, '--k', 'Hz');
h2.DisplayName = 'powV2powI, cap-only';
subplot(2,1,1)
title('Power-amp Voltage to power-amp current');
leg1 = legend([h1, h2]);
set(leg1, 'location', 'northwest')


h3 = frfBode(MF_stage.frf.G_uz2powI, freqs_stage, F2, '-r', 'Hz');
h3.DisplayName = 'uzV2powI, stage';

h4 = frfBode(MF_cap.frf.G_uz2powI, freqs_cap, F2, '--k', 'Hz');
h4.DisplayName = 'uzV2powI, cap-only';

subplot(2,1,1)
title('control Voltage to power-amp current');
leg1 = legend([h3, h4]);
set(leg1, 'location', 'northwest')


