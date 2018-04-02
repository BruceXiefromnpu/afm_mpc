clear
clc

modelFit_file = '/media/labserver/mpc-journal/x-axis_sines_info_out_2-8-2018-01.mat';

load(modelFit_file)
% load /media/labserver/mpc-journal/sysID/x-axis_sines_info_matchTs_out_2-14-2018-01.mat
whos

G_xpow = modelFit.frf.G_xpow_cors;
freqs = modelFit.frf.freq_s;



P_pow2stage = squeeze(G_xpow(1,2, :)./G_xpow(2,2, :));
P_uz2stage = squeeze(G_xpow(1,3, :)./G_xpow(3,3, :));
P_uz2pow = squeeze(G_xpow(2,3, :)./G_xpow(3,3, :));
%%

% Visualize everything
F1 = figure(1); clf
frfBode(P_uz2pow, freqs, F1, 'b', 'Hz');

[P_uz2pow, freqs_pow] = monotonicFRF(P_uz2pow, freqs);
[P_uz2stage, freqs_stage] = monotonicFRF(P_uz2stage, freqs);
[P_pow2stage, freqs] = monotonicFRF(P_pow2stage, freqs);

omegas = freqs*2*pi;
Ts = modelFit.frf.Ts;
P_uz2pow_frd = frd(P_uz2pow, omegas, Ts);

G_uz2pow = tfest(P_uz2pow_frd, 2, 0,5, 'Ts', Ts);
% G_uz2pow.InputDelay = 1
frfBode(G_uz2pow, freqs, F1, '--m', 'Hz');

%%
Vdiv = 2.56/(30.61+2.56);
Vdiv_gain = 1/Vdiv % from resistor measurement
Vdiv_gain = 10.6; % From 9v battery measurement