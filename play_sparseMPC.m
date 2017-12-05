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
x0 = xss*0;

%
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
R = 10000;
if 0
    a = [0 2 2;
        3 3 3;
        4 4 4];
    b = [0 0 1]';
    c = [1 0 0];
    sys = ss(a, b, c, 0, Ts);
    Q1 = eye(3)*2;
    R = 1;
    % Qp = eye(3)*4;
    Qp = dare(a,b,Q1,R);
    mp1 = sparseMPCprob2(sys, 10, Q1, Qp, R);
    x00 = ones(3,1);
else
    S = randn(size(sys_recyc.b,1),1);
   mp1 = sparseMPCprob2(sys_recyc, 100, Q1, Qp, R, S); 
    x00 = SSTools.getXss(sys_recyc);   
   sys = sys_recyc;
end


mp1.add_U_constraint('slew', [-.01, .01])

UX = mp1.solve(x00);

u = UX(1:mp1.nu*mp1.N_mpc);
X = UX(mp1.nu*mp1.N_mpc+1:end);
X = reshape(X, mp1.ns, []);
y_sparse = sys.c*X;


K = dlqr(sys.a, sys.b, Q1, R, S);

syscl = ss(sys.a-sys.b*K, sys.b, sys.c, 0, Ts);

y = initial(syscl, x00);

plot(y)
hold on
plot(y_sparse, '--')







