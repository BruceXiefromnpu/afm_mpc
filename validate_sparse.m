% I want to investigate methods to increase the condition number of the
% Hessian

% Build up the correct model from what is saved from sysID. Ie, put the
% thing in a 
clear
close all
clc
volts2mu = 1;
TOL = 0.01;
trun = 800*40e-6;
ref_f = 1.5;
ref_0 = 0;
umax = 5;

matpath           = getMatPath();
dataroot          = fullfile(matpath, 'AFM_SS', 'System_Identification', 'data','data_xaxis'); 
expName           = ['22-Jun-2016_exp01'];
modFitName    = [expName, '.mat'];
modFitPath    = fullfile(dataroot, modFitName);
load(modFitPath, 'modelFit')


sys = ltiFit(modFitPath, 'SS02').sys;
Nd = 0;
sys.InputDelay = Nd;
Ts = sys.Ts;
PLANT = sys;
% Is this strictly necessary?? Do I need to do it? I think the only
% advantage is that it gives a nice ordering to the modes...
% sysDec    = decimateTF(PLANT_nomH);
% g_drift = sysDec.rz.G(1)*sysDec.rp.G(1)
% 
% sysDec.rz.G(1) = zpk([],[],1, Ts);
% sysDec.rp.G(1) = zpk([],[],1, Ts);
% 


Ns = length(sys.b);

NsNd = Ns+Nd;
sys.InputDelay=Nd;

% Scale system states:
[uss_0, uss_f, ~, ~, xss]   = yss2uss(PLANT, 2, 0);
dcgain_sys = 1/(PLANT.c*xss);
x0 = xss;

%%
% 3). LQR generation gain.        
% -------------------------------------------------------------------------
% -------------------- Constrained LQR Stuff ------------------------------
N_mpc = 12;
slr   = 0.05;

% Pull out open-loop pole-zero information.
[wp_real_x, wz_real_x] = w_zp_real(sys);
rho_1 = wz_real_x(1)/wp_real_x(1);

% zeta_x = [.9, .8, .6, .5 .5];
zeta_x = [.9, .8, .7, .7 .7];
gams_x = [1.5, 1.5, 1.5, 1, 1];
rhos_x = [rho_1*1.0, 1, 1];
 
pint_x = 0.5;
P_x    = getCharDes(sys, gams_x, pint_x, zeta_x, rhos_x, .25);
K_temp = place(sys.a, sys.b, P_x);

%-----------------------------------------------------
[Q0, R1, KKK] = inverseLQR(sys, K_temp);
R1 = 1000;

% Q = blkdiag(Q_nodelay1, zeros(Nd,Nd));

sys_recyc = SSTools.deltaUkSys(sys);
Ns_mpc = size(sys_recyc.B, 1);
[~, Nx, Nu] = getNbar(sys_recyc, [KKK, 0]);
Q1 = blkdiag(Q0, 0);

Qp = dare(sys_recyc.a, sys_recyc.b, Q1, R1); 

mpcProb0 = condensedMPCprob(sys_recyc, N_mpc, Q1, Qp, R1);
mpcProb0.Ainq = [eye(N_mpc); -eye(N_mpc)];
mpcProb0.binq = [zeros(2*N_mpc, 1)+slr];

%%
clc
[Q,R,S,KK] = inverseLQR_cross_weight(sys, K_temp);

%%
clc
clf
% Q1 = eye(3)*2;
% Qp = eye(3)*4;
% R = blkdiag(1);
% a = [0 2 2;
%     3 3 3;
%     4 4 4];
% b = [0 0 1]';
% % 0 2 0]';
% c = [1 0 0];
% sys = ss(a, b, c, 0, Ts);
N = 120;
S = ones(13, 1);
Qp = dare(sys_recyc.a, sys_recyc.b, Q1, R1, S); 
mp1 = sparseMPCprob(sys_recyc, N, Q1, Qp, R1)

UX = mp1.solve([x0;0]);
% UX = mp1.solve([1,1,1]');

% UX = reshape([UX;NaN], 14,[]);

us = UX(end,1:end-1);
Xs = UX(1:end-1, :);
ys = sys_recyc.c*Xs;

plot(ys)

y = lsim(sys_recyc, us, [0:Ts:(length(us)-1)*Ts]', [x0; 0]);
hold on
plot(y, '--')

K = dlqr(sys_recyc.a, sys_recyc.b, Q1, R1, S);
sys_cl = ss(sys_recyc.a-sys_recyc.b*K, sys_recyc.b, sys_recyc.c, 0, Ts)

y2 = initial(sys_cl, [x0;0]);
plot(y2, ':')

%%
% % % % mpcProb1 = FGMprob_fxp(mpcProb0, slr, 20, 32, 30)
% % % mpcProb1 = FGMprob(mpcProb0, slr, 20)
% % % mpcProb1.ziyi = mpcProb1.ziyi + slr;
% % % % mpcProb1.y_i = mpcProb1.z_i;
% % % 
% % % fprintf('condition = %s\n', mpcProb1.kappa)
% % % 
% % % 
% % % mpcParams_fp = struct('K_fdbkLQR', K_temp, 'Q', Q1, 'R', R1, 'TOL', TOL, 'ref_f', ref_f,...
% % %              'ref_0', ref_0, 'N_mpc', N_mpc);
% % % fp_cs = '--k';
% % % xlm = [.004, .033];
% % % ylm = [ref_f - ref_f*TOL*1.1, ref_f + ref_f*TOL*1.1];
% % % 
% % % uss0 = - dcgain_sys*ref_f
% % % z0 = ones(N_mpc,1)*slr;
% % % y0 = z0;
% % % 
% % % sim('LQR_MPC__fictZeros')
% % % 
% % % mpcname = sprintf('MPC: sim (N=%d)', N_mpc);
% % % mpcDist_sim = stepExp_old(y_mpcDist, u_mpcDist, fp_cs, 'MPC: inverse LQR', 'MPC', mpcParams_fp);
% % % F4 = figure(80);
% % % H_mpcDist_sim = plotZoom(mpcDist_sim, F4, xlm, ylm*volts2mu, 'yunits', '[$\mu$m]', 'yscaling', volts2mu, 'textOn', 1);
% % % 
% % % figure(3); 
% % % plot(diff(u_mpcDist.Data))

%%
% Now, compare to the fictitious zeros approach. 

% KKK2 = dlqr(sys.A, sys.B, Q0, 10);
zeros_des = eig(sys.A - sys.B*KKK);
% zeros_des = P_x;
% Move the real, higher frequency pole from the running. Currently, this
% pole is not moved at all from its original location. That means that the
% calculation of (z_i * I - A)^-1 will fail, produce NaN
[~, Ind] = sort(real(zeros_des));
k = find(imag(zeros_des(Ind)) ==0);
zeros_des(Ind(k(1))) = zeros_des(Ind(k(1)))*.95;
% zeros_des(Ind(k(1))) = zeros_des(Ind(k(1)))*.73;
% zeros_des(Ind(k(2))) = zeros_des(Ind(k(2)))*.9997;
zeros_des(Ind(k(2))) = zeros_des(Ind(k(2)))*1.000081% 51;


% [C_hat, D_hat, Q2, S] = getC_hat_S(sys, zeros_des);
C_hat = place(sys.A, sys.B, zeros_des);
D_hat = 1;
Q2 = C_hat'*C_hat;
S = C_hat';
%
sys_temp = ss(sys.a, sys.b, C_hat, D_hat);
fprintf('error in desired vs places polez/zeros %f\n', [sort(zero(sys_temp))- sort(zeros_des)])


R2  = D_hat^2 + 2500;
Q2  = (Q2'+Q2)/2;
Q2 = blkdiag(Q2, 0);
S = [S;0];
Qp2 = dare(sys_recyc.a, sys_recyc.b, Q2, R2); % QQ(end) is s.s. soln to dare.
mpcProb1 = condensedMPCprob_fict(sys_recyc, N_mpc, Q2, Qp2, R2, S);
mpcProb1.Ainq = [eye(N_mpc); -eye(N_mpc)];
mpcProb1.binq = [zeros(2*N_mpc, 1)+slr];

fprintf('condition = %s\n', cond(H3_mpc))


fp_cs = '--k';
xlm = [.004, .033];
ylm = [ref_f - ref_f*TOL*1.1, ref_f + ref_f*TOL*1.1];

uss0 = - dcgain_sys*ref_f
sim('LQR_MPC__fictZeros')

figure(3);clf; plot(diff(u_mpcDist.Data))
mpcname = sprintf('MPC: sim (N=%d)', N_mpc);
mpcDist_sim = stepExp_old(y_mpcDist, u_mpcDist, fp_cs, 'MPC: fict zeros', 'MPC', mpcParams_fp);
F4 = figure(80);
H_mpcDist_sim = plotZoom(mpcDist_sim, F4, xlm, ylm*volts2mu, 'yunits', '[$\mu$m]', 'yscaling', volts2mu, 'textOn', 1);

%%
clc
addpath('~/gradschool/comps/functions')
% LQR Root Locus Plotting.
width = 5;
height = 5;
F5 = figure(5); clf
set(F5, 'color', fig_color,...
    'PaperUnits','inches',...
    'PaperSize',[width height],...
    'PaperPosition',[0,0,width,height],...
    'Units','inches',...
    'Position',[-10,1,width,height]);

[~, hp, hz] = pzplotCL(sys, [], [], F5);

hp.MarkerSize = 12;
hp.LineWidth = 2;
hz.MarkerSize = 12;
hz.LineWidth = 2;

hp = hp(1);
hp.DisplayName = 'Open Loop pole';

hz = hz(1);
hz.DisplayName = 'System zero';

% p1 = plotPdes(zeros_des, Ts, 'g', F5);
p1 = plotPdes(zeros_des, 'color', 'g', 'fig', F5);

p1.Marker = '.';
p1.MarkerSize = 15;
p1.DisplayName = 'fictitious zeros (desired C.L. pole)';

xlim([0.5, 1.1]);
ylim([-.5, .5]);

leg1 = legend([hp, hz, p1]);
leg1.FontSize = 14;
leg1.Position = [0.1432 0.7906 0.7347 0.1729];
set(gca, 'Color', fig_color')
set(gcf, 'InvertHardCopy', 'off');
print(F5, '-dpdf', '~/gradschool/comps/presentation/figures/fict_zero_locus-01.pdf', '-opengl', '-r600')
%%

Q = C_hat'*C_hat;
mpclocus(sys, Q, {D_hat'*D_hat, 1}, 0.01, 10^5, 'nmpc', 24, 'numpoints', 200, 'marker', '.',...
         'color', 'r', 'S', C_hat'*D_hat, 'feedthrough', 1, 'markersize', 15);


xlim([0.5, 1.1]);
ylim([-.5, .5]);


% Plot the placed zeros/ desired poles
% p1 = plotPdes(zeros_des, Ts, 'g', F5)
p1 = plotPdes(zeros_des, 'color', 'g', 'fig', F5)
p1.Marker = '.'
p1.MarkerSize = 15;
%%
% set(gcf, 'InvertHardCopy', 'off');
print(F5, '-dpdf', '~/gradschool/comps/presentation/figures/fict_zero_locus-02.pdf', '-opengl', '-r600')












