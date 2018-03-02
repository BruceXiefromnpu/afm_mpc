clc
Ts = 40e-6;
G = tf([1 -0.55], [1 -1.5 .95], Ts)
pzplot(G)

[G_frf, omegas] = freqresp(G);
G_frf = G_frf(:);

LG = LogCost(Ts, G_frf, omegas);

LG.add_real('z', 1e3)
LG.add_complex('p', .01, 1.2e4)

LG.update_P
LG


G0 = LG.realize();


LG.P
LG.logcost(LG.P)

LG.P = fminsearch(@LG.logcost, LG.P)

% gg = tf(theta(1:2), [1 theta(3:end)], Ts);
gg = LG.realize;
figure(2); clf
bode(G,G0, gg, '--r')


