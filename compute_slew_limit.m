clear



Imax = 100e-3; % 100mA
Ts = 40e-6;
% assert(1/Ts == 25e3);
C = 4.2e-6; % muF

Kamp = 9;  % 20v range to 180 v range

del_Vhigh_max = (Ts/C)*Imax

del_Vlow_max = del_Vhigh_max/Kamp

%%
load /media/labserver/mpc-journal/x-axis_sines_info_out_2-8-2018-01.mat

whos
saveon = 0;
G_xpow = modelFit.frf.G_xpow_cors;
freqs = modelFit.frf.freq_s;



P_pow2stage = squeeze(G_xpow(1,2, :)./G_xpow(2,2, :));
P_uz2stage = squeeze(G_xpow(1,3, :)./G_xpow(3,3, :));
P_uz2pow = squeeze(G_xpow(2,3, :)./G_xpow(3,3, :));


clc
G = ss(modelFit.models.G_uz2pow);
G.InputDelay = 0;
%%
Gc = d2c(G);
Gc = zpk([], pl, 1);
Gc = Gc*(dcgain(G)/dcgain(Gc));

G_vh_I = zpk([0], [], C)*Gc
G_vh_I_z = c2d(G_vh_I, Ts);

step(G_vh_I_z)

bode(G_vh_I_z)
%%
% Gc = d2c(G);
% pl = pole(Gc)
% Gc = zpk([], pl, 1);
% Gc = Gc/dcgain(Gc);
% 
% p1 = -pl(1);
% 
% R1 = 10e3;
% R2 = 5e3;
% C1 = C/150;
% C2 = C/500;
% 
% ff = @(c1c2r1r2)fc1c2(c1c2r1r2, pl);
% 
% c1c2 = fsolve(ff, [C1; C2;R1;R2]);
% 
% C1 = c1c2(1);
% C2 = c1c2(2);
% R1 = c1c2(3);
% R2 = c1c2(4);
% 
% b = 1/C1/R1 + 1/R2/C2 + 1/R2/C1;
% c = 1/(R1*C1*R2*C2); 
% 
% GG = tf(c, [1, b, c])
% pole(GG)
% 
% figure(20)
% bode(Gc, GG)
%%
step(G)
triang = raster(1/(100*Ts), Ts, 400*Ts);
triang.Data = triang.Data*10;
figure(10);
plot(triang.Time, triang.Data)

figure(11);
dtri = diff(triang.Data);
% ./diff(triang.Time);
dtri = [dtri; dtri(end)];
t = triang.Time;
plot(t, dtri)


y = lsim(G, triang.Data, triang.Time);

dy = lsim(G, dtri, t);
figure(12); hold on
plot(t(1:end-1), diff(y))
plot(t, dy, '--')


figure(13); clf
plot(t, dy)
hold on
plot(t, dtri)
xlabel('time [s]')
leg1 = legend('$\Delta y(k)$', '$\Delta u(k)$');
set(leg1, 'FontSize', 14)

save





