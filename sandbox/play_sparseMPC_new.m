% I want to investigate methods to increase the condition number of the
% Hessian

% Build up the correct model from what is saved from sysID. Ie, put the
% thing in a 
clear
close all
addpath('functions')
addpath('models')
clc
volts2mu = 1;
TOL = 0.01;
trun = 800*40e-6;
ref_f = 2;
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


[uss_0, uss_f, ~, ~, xss]   = yss2uss(PLANT, ref_f, 0);
dcgain_sys = 1/(PLANT.c*xss);
x0 = xss*0;

%
% 3). LQR generation gain.        
% -------------------------------------------------------------------------
% -------------------- Constrained LQR Stuff ------------------------------
N_mpc = 8;
du_max   = 0.05;

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
[Q0, R1, K_lqr] = inverseLQR(sys, K_temp);
R1 = 1000;

% Q = blkdiag(Q_nodelay1, zeros(Nd,Nd));

sys_recyc = SSTools.deltaUkSys(sys);
Ns_mpc = size(sys_recyc.B, 1);
[~, Nx, Nu] = getNbar(sys_recyc, [K_lqr, 0]);
Q1 = blkdiag(Q0, 0);
K_lqr = dlqr(sys_recyc.a, sys_recyc.b, Q1, R1);
Qp = dare(sys_recyc.a, sys_recyc.b, Q1, R1); 

mpcProb0 = condensedMPCprob(sys_recyc, N_mpc, Q1, Qp, R1);
mpcProb0.Ainq = [eye(N_mpc); -eye(N_mpc)];
mpcProb0.binq = [zeros(2*N_mpc, 1)+du_max];

%
% clc
% [Q,R,S,KK] = inverseLQR_cross_weight(sys, K_temp);

%%
ref_f = 2.
S = randn(size(sys_recyc.b,1),1);
S = ones(size(sys_recyc.b,1),1);
%
close all
clc
R = 5;
N_mpc = 10;
N_mpc_s = [6, 8, 10, 12];
N_mpc_s = [5,6, 30];
leg_tits = {}
mpc_on=0;
for iter = 1:length(N_mpc_s)
    N_mpc = N_mpc_s(iter);
    
    mp1 = sparseMPCprob(sys_recyc, N_mpc, Q1, Qp, R, S); 
    x00 = SSTools.getXss(sys_recyc);   
    sys = sys_recyc;


    leg_tits{iter} = sprintf('N=%d', N_mpc);
    mp1.add_U_constraint('box', [-du_max, du_max]);

    UX = mp1.solve(x00);

    u = UX(1:mp1.nu*mp1.N_mpc);
    X = UX(mp1.nu*mp1.N_mpc+1:end);
    X = reshape(X, mp1.ns, []);
    y_sparse = sys.c*X;

    if 0
        K = dlqr(sys.a, sys.b, Q1, R, S);
        syscl = ss(sys.a-sys.b*K, sys.b, sys.c, 0, Ts);
        y = initial(syscl, x00);

        plot(y)
        hold on
        plot(y_sparse, '--')
    end
    %
    clc
    trun = 400*Ts;
    mpcProb1 = mp1;
    sim('MPC_fp')
    y1 = y_mpcDist; u1 = u_mpcDist; du1=du_mpcDist;
    

    figure(10); hold on
    plot(y1.Time, y1.Data)
    title('y(t)')

    figure(11); hold on
    plot(u1.Time, u1.Data)
    title('u(t)')
    
    figure(12); hold on
    plot(du1.Time, du1.Data)
    title('du(t)')

    [p1, freqs] = power_spectrum(du1.Data, Ts);
    figure(13); 
    semilogx(freqs, p1, '--')
    hold on
    title('fft(du)')
    
    kcut1 = find(freqs < 400, 1, 'last')
    pwer=[];
    j = 1;
    for k=kcut1:length(freqs)

        p1_wind = p1(k:end);
        L = length(p1_wind);
        pwer(j) = sum(p1_wind.*conj(p1_wind))*L;
        j=j+1;
    end

    figure(14); hold on
    plot(freqs(kcut1:end), pwer, '--')
    drawnow()
    title('bandlimited energy')
    xlabel('freq')
    tilefigs
end

for k=10:14
    
   figure(k)
   legend(leg_tits)
end


%%
mpcProb1 = condensedMPCprob(sys_recyc, N_mpc, Q1, Qp, R, S);

clc
I = eye(N_mpc);
mpcProb1.Ainq = [I; -I];
mpcProb1.binq = [zeros(N_mpc, 1)+du_max];
mpcProb1.binq = [mpcProb1.binq; mpcProb1.binq];
sim('MPC_fp')


figure(10)
hold on
plot(y_mpcDist.Time, y_mpcDist.Data, '--')
title('y(t)')

figure(11)
hold on
plot(u_mpcDist.Time, u_mpcDist.Data, '--')
title('u(t)')

figure(12)
hold on
plot(du_mpcDist.Time, du_mpcDist.Data, '--')
title('delta u(t)')

figure(13)
plot(du1.Data-du_mpcDist.Data)
title('error in delta u(t)')