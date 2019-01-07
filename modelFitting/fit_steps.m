
addpath(fullfile(PATHS.step_exp(), 'functions'))
addpath(fullfile(PATHS.step_exp(), 'functions/canon'))
addpath('hysteresis')
fprintf('============================================\n')
% hyst_file = 'hystID_data_5-4-2018_01.mat';
hyst_file = 'hystID_data_6-1-2018_01.mat';
hyst_path = fullfile(PATHS.sysid, 'hysteresis', hyst_file);
load(hyst_path)

kk = length(hystData.y_exp);
ux = hystData.u_exp(1:kk);
yx = hystData.y_exp(1:kk) - mean(hystData.y_exp(1:500));
tvec = hystData.t_exp(1:kk);

Ts = StageParams.Ts;
du = diff(ux);
idx = find(abs(du) > 0, 1, 'first');

len = idx;
len_fit = 8000;



y = yx(idx-10:end);
u = ux(idx-10:end);
t = tvec(idx-10:end);
figure(10);
N = floor(length(y)/len);
plot(t, u);
hold on
plot(t, y)
%%
Y_mat = [];
U_mat = [];
uss = [];
np = 2;
if np == 1
  zpk0 = [0.9922   0.9997     .8];
  idx_z = [2];
  idx_p = [1];
  idx_k = 3;
%   idx_x0 = [6,7];
elseif np == 2
  zpk0 = [0.9922    0.9997    0.9997    0.9927    .8];
  idx_z = [3,4];
  idx_p = [1,2];
  idx_k = 5;
  % idx_x0 = [6,7];
end
% g0 = ss(zpk(zpk0(idx_z), zpk0(idx_p), zpk0(idx_k), Ts));
% 
% [Nx, Nu] = SSTools.getNxNu(g0)
% x0 = Nx*0;
plants = CanonPlants.plants_ns14();
Gvib = plants.Gvib;
HS = plants.hyst_sat;
% u = PIHyst.sat_op(u, HS.dp, HS.wsp);
% u = PIHyst.hyst_play_op(u, HS.rp, HS.wp, HS.wp*0);
% u = PIHyst.hyst_play_sat_op(u, HS.r, HS.w, HS.d, HS.ws, HS.w*0);

clc
lb = [-1*ones(1, np*2), -Inf]+eps;
ub = -lb;
opts = optimoptions(@lsqnonlin);
opts.MaxFunctionEvaluations = 5000;
opts.MaxIterations = 1000;
opts.Display = 'iter';

theta = [zpk0]; %, x0'];
figure(2);clf
G_s = [];
pls = [];
zrs = [];
for k = 1:N
  
%   x0 = y(1)*Nx;
  theta0 = theta;
%   theta0(idx_x0) = x0;
  
  y_k = y( (k-1)*len+1:k*len);
  u_k = u( (k-1)*len+1:k*len);
  
  Y_mat = [Y_mat, y_k(:)];
  U_mat = [U_mat, u_k(:)];
  

  y_short = y_k(1:len_fit);
  y_short = y_short - y_short(1);
  
  u_short = u_k(1:len_fit);
  u_short = u_short - u_short(1);
  
  t_short = (0:length(u_short)-1)'*Ts;
  u_prime = lsim(Gvib, u_short, t_short);
  
  figure(2); hold on, grid on;
  
  plot(t_short, y_short, 'b')
  plot(t_short, u_short, 'k')
  
  uss(k) = u_k(end);
  Gstatic = zpk([], [], 1, Ts);
%   gdrift_cost = @(theta)fit_gdrift_ss(theta, y_short, u_prime, t_short, idx_p, idx_z, idx_k, idx_x0, Ts);
  gdrift_cost = @(theta)fit_gdrift(theta, Gstatic, y_short, u_prime, t_short, np)
  
%   lb(end-1:end) = 0;
%   ub(end-1:end) = [5, 5]*0;
% theta0 = [0.9922    0.9997    0.9997    0.9927    .8];
  theta = lsqnonlin(gdrift_cost, theta0, lb, ub, opts);
  
  pls(k, :) = theta(idx_p);
  zrs(k, :) = theta(idx_z);
  ks(k)  = theta(idx_k);
  
  g_k = (zpk(theta(idx_z), theta(idx_p), theta(idx_k), Ts));

  y_fitted = lsim(g_k, u_prime, t_short);
  
  plot(t_short, y_fitted, 'r')

  drawnow

  G_s = [G_s, g_k];
  
end
%
%%
figure(4); clf
subplot(4,2,1)
title('pole 1')
hold on, grid on
subplot(4,2,2)
title('pole 2')
hold on, grid on
subplot(4,2,3)
title('zero 1')
hold on, grid on
subplot(4,2,4)
title('zero 2')
hold on, grid on
subplot(4,2,5)
title('gain')
hold on, grid on

dck = ks*0;
for k=1:length(dck)
  dck(k) = dcgain(G_s(k));
end
pls_w = (-log(pls)/Ts)/2/pi;
zrs_w = (-log(zrs)/Ts)/2/pi;

subplot(4,2,1)
plot(uss, pls_w(:, 1), 'x')
plot(uss, pls_w(:,1))
ylabel('Hz')
xlabel('$u_{ss}$')

subplot(4,2,3)
plot(uss, zrs_w(:, 1), 'x')
plot(uss, zrs_w(:,1))
ylabel('Hz')
xlabel('$u_{ss}$')

subplot(4,2,5)
plot(uss, dck, 'x')
plot(uss, dck)
xlabel('$u_{ss}$')


if np ==2
  subplot(4,2,2)
  plot(uss, pls_w(:, 2), 'x')
  plot(uss, pls_w(:,2))
  ylabel('Hz')
  xlabel('$u_{ss}$')

  % plot(uss, pls1 - zrs1)

  subplot(4,2,4)
  plot(uss, zrs_w(:, 2), 'x')
  plot(uss, zrs_w(:,2))
  ylabel('Hz')
  xlabel('$u_{ss}$')
  % plot(uss, (pls1+zrs1)/2)

    subplot(4,2,6)
  plot(uss, zrs_w(:,1) - pls_w(:,2), 'b')
  hold on, grid on
  plot(uss, zrs_w(:,1) - pls_w(:,2), 'x')
  title('z1 - p2')

  subplot(4,2,7)
  plot(uss, (zrs_w(:,2) - pls_w(:,1)), 'b')
  hold on, grid on
  plot(uss, (zrs_w(:,2) - pls_w(:,1)), 'x')
  title('(z2 - p1)')
  ylabel('Hz')
  
  subplot(4,2,8)
  plot(uss, (pls_w(:,1) + zrs_w(:,2))/2, 'b')
  hold on, grid on
  plot(uss, (pls_w(:,1) + zrs_w(:,2))/2, 'x')
  title('(p1 + z2)/2')

else
    subplot(4,2,7)
    plot(uss, (zrs_w(:,1) -pls_w(:,1) ), 'b')
    hold on, grid on
    plot(uss, (zrs_w(:,1) - pls_w(:,1) ), 'x')
    title('(z1 - p1)')
    ylabel('Hz')

end


%%
save('../many_steps_fits.mat', 'G_s', 'pls', 'zrs', 'ks')
















