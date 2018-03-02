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


% Visualize everything
F1 = figure(1); clf
frfBode(P_pow2stage, freqs, F1, 'g', 'Hz');
frfBode(P_uz2stage, freqs, F1, 'r', 'Hz');
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
% Try fitting everything at once, (uz 2 stage).
G = P_uz2stage;

F2 = figure(2); clf;
frfBode(P_uz2stage, freqs, F2, 'r', 'Hz');
% frfBode(G, freqs, F, 'g', 'Hz')
ss_opts = frf2ss_opts('Ts', Ts);

f2ss = frf2ss(G, omegas, 11, ss_opts); % 12


sys = f2ss.realize(12); % 12
nmp_z = tzero(sys);
nmp_z = nmp_z(find(abs(nmp_z) > 1));
g_eject = zpk([], nmp_z(1:2), 1, Ts);
sys = minreal(g_eject*sys)/dcgain(g_eject);
sys.InputDelay = 9;
gr = zpk([], exp(-2*pi*2500*Ts), 1, Ts);
gr = gr/dcgain(gr);


frfBode(sys*gr, freqs, F2, '--k', 'Hz');

plotPZ_freqs(sys*gr, F2);
%
% grab g_drift;
pl = pole(sys);
zr = zero(sys);
pl_real = pl(find(imag(pl)==0));
zr_real = zr(find(imag(zr)==0));
gdrift = zpk(zr_real(1), pl_real(end), 1, Ts);
% gdrift = gdrift/dcgain(gdrift);

%%
clc
nd = 11;
kk = 65;
ll = 20;
Gdelay = zpk([], zeros(1, nd), 1, Ts);
Gdelay_frf = freqresp(Gdelay, omegas);
P_uz2stage_nod = P_uz2stage./Gdelay_frf;
P_lowfreq_frd = frd(P_uz2stage_nod(ll:kk), omegas(ll:kk), Ts);
G_lowfreq = tfest(P_lowfreq_frd, 2, 2, 'Ts', Ts, 'FeedThrough', 1);
zpk(ss(G_lowfreq))

figure(20); clf
H = bodeplot(gdrift, G_lowfreq,'-.', frd(P_uz2stage, omegas), 'r', P_lowfreq_frd, '--', g)
setoptions(H, 'FreqUnits', 'Hz')
grid
%%
hold on
bodeplot()

%%

% Try fitting the pow2stage and uz2pow individually, and multiplying
% together.

F3 = figure(3);clf
frfBode(P_pow2stage, freqs, F3, 'r', 'Hz');
frfBode(P_uz2stage, freqs, F3, 'g', 'Hz');
frfBode(P_uz2pow, freqs, F3, 'b', 'Hz');
frfBode(G_uz2pow, freqs, F3, '--m', 'Hz');


G = P_uz2stage./squeeze(freqresp(G_uz2pow, omegas));
% G = P_pow2stage;
% frfBode(G, freqs, F, 'g', 'Hz')
ss_opts = frf2ss_opts('Ts', Ts)

f2ss = frf2ss(G, omegas, 6, ss_opts);


sys = f2ss.realize(11); %% 12
nmp_z = tzero(sys);
nmp_z = nmp_z(find(abs(nmp_z) > 1))
% g_eject = zpk([], nmp_z(1:2), 1, Ts)
% sys = minreal(g_eject*sys)/dcgain(g_eject)
frfBode(sys, freqs, F3, '--k', 'Hz');



plotPZ_freqs(sys, F3)

F4 = figure(4);clf
% frfBode(P_pow2stage, freqs, F3, 'r', 'Hz');
frfBode(P_uz2stage, freqs, F4, 'g', 'Hz');
% frfBode(P_uz2pow, freqs, F3, 'b', 'Hz');
% frfBode(G_uz2pow, freqs, F3, '--m', 'Hz');
gdrift = 1;
frfBode(sys*G_uz2pow*gdrift, freqs, F4, '--k', 'Hz');
plotPZ_freqs(sys*G_uz2pow*gdrift, F4);

%%
frfBode(sys*G_uz2pow, freqs, F1, 'r', 'Hz')

plotPZ_freqs(sys*G_uz2pow, F1)

%%

% modelFit.models.G_uz2pow = G_uz2pow;
% modelFit.models.G_pow2uz = sys;

modelFit.models.G_uz2stage = sys;
modelFit.models.G_real_extra = gr;
save(modelFit_file, 'modelFit')
% save('/media/labserver/mpc-journal/sysID/x-axis_sines_info_matchTs_out_2-14-2018-01.mat', 'modelFit')
%%










