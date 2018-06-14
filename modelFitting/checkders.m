
ts = 0.01;
t = (0:1000)'*ts
u = sin(2*pi*t);
gg = c2d(tf(100, [1, 20, 100]), ts);
y = lsim(gg, u, t)

Jfunc = @(theta) hyst_drift_paralel_der(theta, u, y, t, Nhyst, Nsat, nsd, gg, r, d);

clc
% u = sin(2*pi(0:0.01
dl =0.00001;
idx = 16;
DL = theta0*0;
DL(idx) = DL(idx)+dl;

[err, udrift1, Derr] = Jfunc(theta0);

[err_plus_dl, udrift2] = Jfunc(theta0+DL);

der_fin = (err_plus_dl - err)/dl;



figure(10); clf
subplot(2,2,1)
plot(der_fin)
title('finite diff')

subplot(2,2,2)
plot(Derr(:,idx))
title('provided der')

subplot(2,2,3)
plot(Derr(:, idx) - der_fin)
title('der error')

subplot(2,2,4)
plot(der_fin)
hold on
plot(Derr(:,idx), '--')
title('both ders')

figure(11); clf
plot(udrift1)
hold on
plot(udrift2, '--')