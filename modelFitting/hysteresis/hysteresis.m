clc

% %%
% clc
Ts = 0.1;
t1 = 15;
tf = 100;
umax = 3;
k1 = 0.09;

u_reset = PIHyst.gen_reset_u(t1, tf, Ts, k1, umax);

% figure(5)
% plot(T, u_reset)

u0 = (0:0.1:3-0.1);
u = (3-0.1:-0.1:-3);
u = [u, fliplr(u)];
u = [u0, repmat(u, 1, 3), fliplr(u0)];
u = [u'; u_reset];
% u = diff(u);
% u = [u, u(end)]*10;

r = [0, .5, 1, 1.5, 2, 2.5, 3]*2;

% w = [.1,.5, 1, 1.25, 1.5, 1.6, 1.75]';
w = [1,1, 1, 1, 1, 1, 1]';
% w = w/sum(w);
% y0 = repmat(0, length(r), 1);

[y, y_vec] = hyst_play_op(u, r, [w; y0]);
fprintf('yfinal = %.5f\n', y(end));
fprintf('y_vec final = %5f\n', y_vec(end, :))
figure(1); clf;
plot_ind(u', y', '-or')
grid on
hold on;
plot(u(end), y(end), 'x')
xlm = xlim;

figure(2)
plot(u)
grid on


%%
fname = fullfile(PATHS.sysid, 'hyst_id_4-20-2018_01.csv');

dat = csvread(fname);
kstart = 1;
u = dat(kstart:end, 1);
y = dat(kstart:end, 2);
% y = y - mean(y);
y = y - y(1);
figure(10)
plot(u)
%%


y0 = repmat(0, length(r), 1);
w_y0_0 = [w; y0];

fun = @(wy0) hyst_play_op(u, r, wy0) - y;
w_y0 = lsqnonlin(fun, w_y0_0);

[w_y0_0, w_y0]

y_est = hyst_play_op(u, r, w_y0);

figure(2); clf
plot(u, y)
grid on, hold on;
plot(u, y_est, '--r')

figure(20)
plot([y, y_est])


%%
fname = fullfile(PATHS.sysid, 'hyst_id_4-20-2018_02.csv');

dat = csvread(fname);
kstart = 1;
u2 = dat(kstart:end, 1);
y2 = dat(kstart:end, 2);
V2 = dat(kstart:end, 3);
I2 = dat(kstart:end, 4)/15;
y2 = y2 - y2(1);
%%
clc
n = length(r);
w = w_y0(1:n);
% fun2 = @(y0) hyst_play_op(u2, r, [w; y0]) - y2;
fun2 = @(wy0) hyst_play_op(u2, r, wy0) - y2;
wy02 = lsqnonlin(fun2, w_y0);


%%
% wy03 = [w_y0(1:n); wy02(n+1:end)];
figure(3); clf
plot(y2)
hold on, grid on;
plot(u2*0.725, '-k', 'LineWidth', 2)
y2_est = hyst_play_op(u2, r, wy02);

plot(y2_est, '--r')


%%

fun3 = @(w_y01_y02) [hyst_play_op(u, r, w_y01_y02(1:2*n)) - y;
                    hyst_play_op(u2, r, w_y01_y02([1:n, 2*n+1:3*n]')) - y2];                    

w_y01_y02 = lsqnonlin(fun3, [w_y0; wy02(n+1:end)]);
%%
w_y01 = w_y01_y02(1:2*n);
w_y02 = w_y01_y02([1:n, 2*n+1:3*n]')
w_y02(n+1:end) = 0;

y1_est = hyst_play_op(u, r, w_y01);
y2_est = hyst_play_op(u2, r, w_y02);

figure(100); clf
subplot(3,1,1)
plot(u, y)
hold on, grid on
plot(u, y1_est)

subplot(3,1,2)
plot(u2, y2)
hold on, grid on
plot(u2,y2_est)
plot(u2, u2)

subplot(3,1,3)
plot(y2)
hold on, grid on
plot(y2_est)
%%

fun3 = @(w) [hyst_play_op(u, r, w) - y; hyst_play_op(u2, r, w) - y2];                    

w_y01_y02 = lsqnonlin(fun3, w_y0(1:n));
%%
clc
w_y01 = [w_y01_y02; zeros(n, 1)];
w_y02 = [w_y01_y02; zeros(n, 1)];
% w_y02(n+1:end) = 0;

y1_est = hyst_play_op(u, r, w_y01);
y2_est = hyst_play_op(u2, r, w_y02);
%%
clc
[rp, wp] = invert_hyst_PI(r, w_y02(1:n));

u2_inv = hyst_play_op(u2, rp, [wp; zeros(n, 1)]);
u2_prime = hyst_play_op(u2_inv, r, w_y02);
figure(300); clf
plot(u2), hold on, plot(u2_prime, '--')
plot(u2_inv, ':')
%%


figure(100); clf
% subplot(2,1,1)
plot(u, y)
hold on, grid on
plot(u, y1_est)


% subplot(2,1,2)
figure(200); clf
plot(y2)
hold on, grid on
plot(y2_est)
plot(u2*.68)
plot(u2_inv*.68)
%%
clc
Ts = 0.01;
t1 = 5;
tf = 10;
T1 = (1:t1/Ts)'*Ts;
T2 = (t1/Ts+1:tf/Ts)'*Ts;
T = [T1; T2];

umax = 10;
k1 = 1;
k2 = 1;

a1 = umax - k1*T1;
a2 = a1(end)*exp(k2*(t1-T2));
a = [a1; a2];

omega = 5*2*pi;
u = a.*sin(omega*T - pi/2);

figure(5)
plot(T, u)











