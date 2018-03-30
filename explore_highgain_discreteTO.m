% ----------------------------------------------------------------------- %
% ----------------------- Create a toy model ---------------------------- %
clc
Ts = 5e-6;

w1 = 650*2*pi;
z1 = 0.02;
w2r = 500*2*pi;
s1 = -w1*z1 + 1j*w1*sqrt(1-z1^2);
s2 = conj(s1);
s3 = -w2r;
% real poles
ps = [s1, s2, s3];
% discrete poles
pz = exp(ps*Ts);
G = zpk([], pz, 1, Ts);
sys = ss(G)/dcgain(G);

[Gfrf, ws] = freqresp(sys, logspace(0, 3, 100));
Gfrf = Gfrf(:);
F1 = figure(20);
frfBode(sys, ws*2*pi, F1, '-r', 'Hz');
plotPZ_freqs(sys, F1);


sys_nodelay = sys;
PLANT = sys;
Ns = length(sys.b);
Nd = 0;

xss = SSTools.getNxNu(PLANT);
dcgain_sys = 1/(PLANT.c*xss);
x0 = xss*0;


% 3). LQR generation gain.        
% -------------------------------------------------------------------------
% -------------------- Constrained LQR Stuff ------------------------------
sys_recyc = SSTools.deltaUkSys(sys);
Ns_mpc = size(sys_recyc.B, 1);
[Nx, Nu] = SSTools.getNxNu(sys_recyc);


zeta_x = [.95];
gams_x = [2.5];
rhos_x = 2.5;
 

% P_x    = getCharDes(sys_nodelay, gams_x, [], zeta_x, rhos_x, .25);
% K_temp = place(sys_nodelay.a, sys_nodelay.b, P_x);
% [Q0, R1, K_lqr] = inverseLQR(sys_nodelay, K_temp);
P_x    = getCharDes(sys_recyc, gams_x, [.7], zeta_x, rhos_x, .25);
P_x = [0.6729 + 0.0863i   0.6729 - 0.0863i,...
      0.625 + 0.2j, 0.625-0.2j];
K_temp = place(sys_recyc.a, sys_recyc.b, P_x);
[~,~,Q0, S, Rmin] = fict_zeros(sys_recyc, P_x);

Q1 = Q0;
S1 = S;

% K_lqr = dlqr(sys_recyc.a, sys_recyc.b, Q1, 357, S1);


%%

sys_sim = sys_recyc;
[Nx_sim, Nu_sim] = SSTools.getNxNu(sys_sim);
x0_sim = Nx_sim*0;


du_max = 0.5;
xf = Nx_sim*5.5;
k0 = 93;

toBisect = TimeOptBisect(sys_sim, du_max);
toBisect.max_iter = 50;
toBisect.TOL = 1e-6;


[X, U, status]=toBisect.time_opt_bisect(x0_sim, xf, 'k0', k0);

Nopt = length(U.Data)

u = [U.Data; ones(2*k0,1)*ref_k*Nu_sim]; % zero for deltaUk
t = [0:1:length(u)-1]'*sys_sim.Ts;

[y, t, x] = lsim(sys_sim, u, t, x0_sim);



%
% clc
figure(1);clf
subplot(2,1,1)
plot(y)
hold on, grid on

subplot(2,1,2)
plot_ind(u)
hold on, grid on

%%
N = length(U.Data);

xerr = -X.Data +repmat(xf', N, 1);
% xerr = -x +repmat(xf', size(x,1), 1);

k_satp = find(u >= 0.5-.0001)
k_satm = find(u <= -0.5 +.0001)

% k_usat = [11, 18, 29];
uopt = u(1:Nopt);
k_usat = find(uopt <=0.5-0.0001 & uopt>=-0.5+0.0001);
% k_usat = 1:Nopt+10;

% x_usat = xf - X.Data(k_usat', :)
x_usat = xerr(k_usat', :)
u_usat  = U.Data(k_usat)
% u_usat = u(k_usat);

X_satp = xerr(k_satp, :);
X_satm = xerr(k_satm, :);
u_satp = u(k_satp);
u_satm = u(k_satm);

X_sat = [-X_satp;
         X_satm];
u_sat = [-repmat(0.5, length(k_satp), 1);
          repmat(-0.5, length(k_satm), 1)];
      

% K_est = (x_usat\u_usat)'
K_est = lsqlin(x_usat, u_usat, X_sat, u_sat)';
% K_est = lsqlin(x_usat, u_usat)';




uvec_est = zeros(N, 1);
yvec_est = zeros(N, 1);
xvec_est = zeros(4, N);
xvec_est(:, 1) = xerr(1,:);


uvec_lqr = zeros(N, 1);
yvec_lqr = zeros(N, 1);
xvec_lqr = zeros(4, N);
xvec_lqr(:, 1) = xerr(1,:);

%
xk = xerr(1,:)'*0;
xk_lqr = xerr(1,:)'*0;
for kk = 1:N*2
%     yvec(kk) = sys_sim.c*xk;
    xvec_est(:,kk) = xk;
    uk = K_est*(xf - xk);
    uk = max(min(uk, du_max), -du_max);
    uvec_est(kk) = uk;
    xk = sys_sim.a*xk + sys_sim.b*uk;

    % -------------- Now LQR ---------------%
    xvec_lqr(:,kk) = xk_lqr;
    uk_lqr = K_lqr*(xf - xk_lqr);
    % (norm(K_est)/norm(K_lqr))*
    uk_lqr = max(min(uk_lqr, du_max), -du_max);
    uvec_lqr(kk) = uk_lqr;
    xk_lqr = sys_sim.a*xk_lqr + sys_sim.b*uk_lqr;

    
    
end


figure(1);
subplot(2,1,1)
hold on
plot(sys_sim.c*xvec_est, '--')
% plot(sys_sim.c*xvec_lqr, ':')


subplot(2,1,2)
hold on
uvec_est = max(min(uvec_est, du_max),-du_max);
plot(uvec_est, '--')
plot(k_usat, u(k_usat), 'x')
% plot(uvec_lqr, ':')



% plot(k_usat, u_sat, 'xk')
%%
kkdot = -K_lqr*K_est'
costheta = norm(K_lqr, 2)*norm(K_est, 2)/kkdot
acosd(costheta)



