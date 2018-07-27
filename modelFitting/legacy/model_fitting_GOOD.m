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
frfBode(P_pow2stage, freqs, F1, 'Hz', 'g');
frfBode(P_uz2stage, freqs, F1, 'Hz', 'r');
frfBode(P_uz2pow, freqs, F1, 'Hz', 'b');

[P_uz2pow, freqs_pow] = monotonicFRF(P_uz2pow, freqs);
[P_uz2stage, freqs_stage] = monotonicFRF(P_uz2stage, freqs);
[P_pow2stage, freqs] = monotonicFRF(P_pow2stage, freqs);

omegas = freqs*2*pi;
Ts = modelFit.frf.Ts;
P_uz2pow_frd = frd(P_uz2pow, omegas, Ts);

G_uz2pow = tfest(P_uz2pow_frd, 2, 0,5, 'Ts', Ts);
% G_uz2pow.InputDelay = 1
frfBode(G_uz2pow, freqs, F1, 'Hz', '--m');

%%
% Try fitting everything at once, (uz 2 stage).
G = P_uz2stage;

F2 = figure(2); clf;
frfBode(P_uz2stage, freqs, F2, 'Hz', 'r');
% frfBode(G, freqs, F, 'Hz', 'g')
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


frfBode(sys*gr, freqs, F2, 'Hz', '--k');

plotPZ_freqs(sys*gr, F2);
%
% grab g_drift;
pl = pole(sys);
zr = zero(sys);
pl_real = pl(find(imag(pl)==0));
zr_real = zr(find(imag(zr)==0));
%%
gdrift = zpk(zr_real(1), pl_real(end), 1, Ts);
gdrift = gdrift/dcgain(gdrift);



% Try fitting the pow2stage and uz2pow individually, and multiplying
% together.

F3 = figure(3);clf
frfBode(P_pow2stage, freqs, F3, 'Hz', 'r');
frfBode(P_uz2stage, freqs, F3, 'Hz', 'g');
frfBode(P_uz2pow, freqs, F3, 'Hz', 'b');
frfBode(G_uz2pow, freqs, F3, 'Hz', '--m');


G = P_uz2stage./squeeze(freqresp(G_uz2pow, omegas));
% G = P_pow2stage;
% frfBode(G, freqs, F, 'Hz', 'g')
ss_opts = frf2ss_opts('Ts', Ts)

f2ss = frf2ss(G, omegas, 6, ss_opts);


sys = f2ss.realize(11); %% 12
nmp_z = tzero(sys);
nmp_z = nmp_z(find(abs(nmp_z) > 1))
% g_eject = zpk([], nmp_z(1:2), 1, Ts)
% sys = minreal(g_eject*sys)/dcgain(g_eject)
frfBode(sys, freqs, F3, 'Hz', '--k');



plotPZ_freqs(sys, F3)

F4 = figure(4);clf
% frfBode(P_pow2stage, freqs, F3, 'Hz', 'r');
frfBode(P_uz2stage, freqs, F4, 'Hz', 'g');
% frfBode(P_uz2pow, freqs, F3, 'Hz', 'b');
% frfBode(G_uz2pow, freqs, F3, 'Hz', '--m');
gdrift = 1;
frfBode(sys*G_uz2pow*gdrift, freqs, F4, 'Hz', '--k');
plotPZ_freqs(sys*G_uz2pow*gdrift, F4);

%%
frfBode(sys*G_uz2pow, freqs, F1, 'Hz', 'r')

plotPZ_freqs(sys*G_uz2pow, F1)

%%
clc
w_s = omegas;
Nd = 6;
gdelay = (exp(1i*w_s*Ts)).^(Nd+0);
Gfrf = G.*gdelay;


% Set number freqs and q, 

q = 80;
p = 1;  % no ouputs ?
n = 12;


opts = struct('n', n, 'q', q, 'p', p, 'Feedthrough', 0);
R_s = ones(length(w_s), 1)*1;
R_s = linspace(.001, 1, length(w_s) );
% R_s = logspace(log10(0.01), 0, length(w_s));


colrs = {'k', 'g', 'm', 'b', 'y'};
% for i = 1:10
% k = 280;
% [A, B, C, D] = frf_id_nonunif_mckelvey(Gfrf(1:k), w_s(1:k)*Ts, R_s(1:k), opts);

[A, B, C, D] = frf_id_nonunif_mckelvey(Gfrf, w_s*Ts, R_s, opts);

sys2 = ss(A, B, C, D*0, Ts, 'InputDelay', Nd);
nmp_z = tzero(sys2);
nmp_z = nmp_z(find(abs(nmp_z) > 1))
% g_eject = zpk([], nmp_z(1:2), 1, Ts)
% sys2 = minreal(g_eject*sys2)/dcgain(g_eject)



F4 = figure(4); clf
frfBode(G, freqs, F4, 'Hz', 'r')
frfBode(sys2, freqs, F4, 'Hz', '--g');

% frfBode(P_uz2stage, freqs, F2, 'Hz', 'k');
plotPZ_freqs(sys2, F4);
% frfBode(sys2*G_uz2pow, freqs, F2, 'Hz', '--g');
% plotPZ_freqs(sys2*G_uz2pow, F2);



%%
if 0
    F3 = figure(3); clf;
    frfBode(P_uz2stage,  freqs, F3, 'Hz', 'r');
    frfBode(sys2*G_uz2pow, freqs, F3, 'Hz', '--k');

    plotPZ_freqs(sys2*G_uz2pow, F3);


    sys_fit_pow_through_stage = sys2*G_uz2pow;
    modelFit.models.G_uz2pow = G_uz2pow;
    modelFit.models.G_pow2stage = sys2;
end
modelFit.models.G_uz2stage = sys2;


%%
% modelFit.models.G_uz2pow = G_uz2pow;
% modelFit.models.G_pow2uz = sys;

modelFit.models.G_uz2stage = sys;
modelFit.models.G_real_extra = gr;
save(modelFit_file, 'modelFit')
% save('/media/labserver/mpc-journal/sysID/x-axis_sines_info_matchTs_out_2-14-2018-01.mat', 'modelFit')
%%










