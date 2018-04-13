clc, clear
C_stage = 3.8e-6;
Imax = 100e-3;
Ts = 40e-6;


Vdiv = 2.56/(19.61+2.56);
% Vdiv_gain = 1/Vdiv; % from resistor measurement
Vdiv_gain = 10.6; % From 9v battery measurement
load( fullfile(PATHS.exp, 'x-axis_sines_info_out_2-8-2018-01.mat'))
Gpow_hat = modelFit.models.G_uz2pow;

Gpow = ss(Gpow_hat*Vdiv_gain);
Gpow.InputDelay = 0;
Gpow = minreal(Gpow*zpk([], [0 0], 1, Ts))


dVmax = (Ts/C_stage)*Imax
dVmax = 0.5;

G_stage = modelFit.models.G_pow2uz/Vdiv_gain;

G_stage.InputDelay = 0;
sys_old = modelFit.models.G_uz2stage;
sys = G_stage*Gpow;
% sys = modelFit.models.G_uz2stage;
F1 = figure(1000);

[~,~,omegas] = bode(sys);
frfBode(Gpow, omegas/2/pi, F1,'-b', 'Hz');
frfBode(G_stage, omegas/2/pi, F1,'-k', 'Hz');
frfBode(sys_old, omegas/2/pi, F1,'-r', 'Hz');
frfBode(sys, omegas/2/pi, F1,'-g', 'Hz');

plotPZ_freqs(sys, F1);
plotPZ_freqs(G_stage, F1);

% figure(100)
% pzplot(Gpow, 'b', G_stage, 'k')
% 
% figure(3)
% pzplot(sys_old, 'g', sys, 'r')
% sys_nodelay = sys;
sys_nodelay = modelFit.models.G_uz2stage;
sys2 = modelFit.models.G_uz2stage;
[~, gdrift2] = eject_gdrift(sys2);

sys_nodelay.InputDelay = 0;
% [sys_nodelay, gdrift] = eject_gdrift(sys_nodelay)

%%
% Now, try the state constraint with the condensed formulation.
x0_pow = [0;0];
du_max = 2.0;

ref = 1; %OA works, QP fails
Nx = SSTools.getNxNu(sys_nodelay);
x0 = -Nx*ref;


% ----------------------------------------------------------------------- %
% -------------------- Now, try the time-optimal ------------------------ %

clc
% k0 = 108 % ref=4
k0 =78


CON3 = CondenCon(Gpow, x0_pow, k0);
dVmax2 = 0.3
du_max2 = 2;
CON3.add_state_con('box', dVmax2*0.8);
CON3.add_input_con('box', du_max2);
CON3.update_binq

[sys_TO, gdrift] = eject_gdrift(sys_nodelay);
% sys_TO = sys_nodelay;
sys_TO = SSTools.deltaUkSys(sys_TO);
Nx_TO = SSTools.getNxNu(sys_TO);
%
x0_TO = Nx_TO*0;
xf = Nx_TO*ref*1;
P = path;
addpath(genpath(fullfile(getMatPath, 'solvers/cvx')))
C=[];
%
% du_max = 0.2;
for k=0:k0-1
    C = [sys_TO.a^k*sys_TO.b C];
end
e_k0 = zeros(1,k0);
e_k0(end) = 1; % Kth unit vector. [0 0 .... 0 0 1]. To pick out last element of u.
I = eye(length(sys_TO.b));

cvx_begin
    variable u(k0);
%     variable t;
    L = sys_TO.a^k0*x0_TO + C*u;
    minimize norm(L - xf)
    % minimize norm(u)
    % minimize norm(t)
    subject to
%     norm(u, Inf) <= du_max;
    [CON3.Ainq; -CON3.Ainq]*u <= [CON3.ubAinq; -CON3.lbAinq];
    u <= CON3.ub;
    -u<= -CON3.lb;
cvx_end
path(P);
% rmpath(genpath(fullfile(getMatPath, 'solvers/cvx')))

fprintf('Results optval = %f\n', cvx_optval)
results.U = u;
results.cvx_optval = cvx_optval;
results

u = [u; repmat(0, 400, 1)];

t = [0:1:length(u)-1]'*Ts;
utraj_inv = lsim(gdrift, u, t);
uref = cumsum(utraj_inv);
% uref = cumsum(u);
%
sys_delay = sys_nodelay;
y = lsim(sys_nodelay, uref, t, sys_nodelay.b*0);
sys_delay.InputDelay = 10;
ydelay = lsim(sys_delay, uref, t, sys_delay.b*0);
ypow = lsim(Gpow, utraj_inv, t, x0_pow*0);
%
F1 = figure(9); clf;
subplot(3,1,1), hold on
h5 = plot(t, y, '-k');
h5.DisplayName = 'Time-Optimal'
grid on

% legend([h3,h4,h5]);
subplot(3,1,2), hold on
plot(t, u, '-k')
plot(t(k0), u(k0), 'x')
title('$\Delta u$')
grid on


subplot(3,1,3), hold on, grid on
plot(t, cumsum(u), '-k')

title('u(k)')

figure(200), clf
hold on
h6 = plot(ypow);
h6.DisplayName = 'TO sim';

% legend([h1, h2, h6])
xlim([0, 300])
    grid on
    title('Power Amplifer $\Delta y$')
    xlm = xlim;
    hold on;
    plot(xlm, [dVmax, dVmax], '--k')
    plot(xlm, -[dVmax, dVmax], '--k')



%%
% ---- Paths for shuffling data to labview and back. ------
controlParamName = 'exp01Controls.csv'; 
refTrajName      = 'ref_traj_track.csv';
outputDataName = 'exp01outputBOTH.csv';
%labview reads data here
controlDataPath = fullfile(PATHS.step_exp, controlParamName); 
% labview saves experimental results/data here
dataOut_path    = fullfile(PATHS.step_exp, outputDataName); 
% labview reads desired trajectory here
refTrajPath     = fullfile(PATHS.step_exp, refTrajName); 
% location of the vi which runs the experiment.
% vipath ='C:\Users\arnold\Documents\labview\ACC2017_archive\play_AFMss_integral_trajTrack.vi';
% vipath ='C:\Users\arnold\Documents\MATLAB\afm_mpc_journal\labview\play_AFMss_integral_trajTrack_nod.vi';
vipath ='C:\Users\arnold\Documents\MATLAB\afm_mpc_journal\labview\play_AFM_PI_trajTrack.vi';


%----------------------------------------------------
% Save the controller to .csv file for implementation
clc
vipath ='C:\Users\arnold\Documents\MATLAB\afm_mpc_journal\labview\play_AFM_PI_trajTrack.vi';
clear vi; clear e;
% delay before we start tracking, to let any transients out. Somewhere of a
% debugging setting. 
SettleTicks = 20;  
Iters  = 400
Iters = min(Iters, length(uref)-10);


% creat and pack data. Then save it. 
tt = t;
% ref_traj.Time;
yy = ydelay;
% ref_traj.Data;
% uKx  = yy*Nbar;
Nbar = 1;
[~, ~, y_uKx] = pack_uKx_y(uref, yy, tt);

% AllMatrix = packMatrix(sys_obs, L, K);
saveControlDataPITrack(y_uKx, refTrajPath);
%
% -----------------------RUN THE Experiment--------------------------------
[e, vi] = setupVI(vipath, 'SettleTicks', SettleTicks, 'Iters', Iters, 'Stop', 0,...
            'Ki', 0.02*0, 'umax', 5, 'outputData BOTH', dataOut_path,...
            'reference traj path', refTrajPath);
vi.Run
% -------------------------------------------------------------------------

        
% Now, read in data, and save to structure, and plot.
AFMdata = csvread(dataOut_path);
%
clc
[y_exp, u_exp, ypow_exp] = unpackExpDataPITrack(AFMdata, Ts);
size(y_exp.Time)
figure(10), clf
subplot(2,1,1)
hold on
plot(t, ydelay)
plot(y_exp.Time, y_exp.Data, '--r')
grid on

subplot(2,1,2)
hold on
plot(t, uref)
plot(u_exp.Time, u_exp.Data, '--r')
grid on
legend('sim', 'exp')

figure(200)
dypow = diff(ypow_exp.Data)*Vdiv_gain;
h7 = plot(dypow(3:end), '--r');
h7.DisplayName = 'exp';

legend([h6, h7])

figure(400); clf
plot(diff(ypow))
hold on
plot(diff(dypow(3:end)))
title('$\Delta^2 y_{pow}$')
%%
sysinv = 1/zpk(sys_nodelay);
sysinv.InputDelay = 3;
sysinv = absorbDelay(sysinv)
ye = lsim(sysinv, y, t);
figure, plot(t, ye, t, y, '--')

%%


expOpts = stepExpOpts(linOpts, 'pstyle', 'k', 'name',  'AFM Stage');

afm_exp = stepExp(y_exp, u_exp, expOpts);
H2 = plot(afm_exp, F1);



