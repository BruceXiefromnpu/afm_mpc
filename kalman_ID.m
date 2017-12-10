clc
% close all
clear
% load('/home/arnold/matlab/publications/afmMPC/feb27_2017_exp01_yr_1p0.mat')
addpath('functions')
addpath(genpath('../vanoverschee'))
dat = csvread('/media/labserver/sysID/out_12_06.csv');
dat = dat';
size(dat)

dataPath = fullfile(getMatPath, 'publications', 'afmMPC');

modelFitPath = '/media/labserver/sysID/x-axis_sines_info_out_12-10-2017-04.mat';
load(modelFitPath);
w_s_x = modelFit.frf.w_s;
Gx_frf = modelFit.frf.G_frf;
% load(fullfile(dataPath, 'mimo_modelData.mat'))
%%
u = dat(2:end, 1);
y = dat(2:end, 2);
y = y - mean(y);
u = u - mean(u);
y = y(1000:end);
u = u(1000:end);
% plot(dat(:,1))



% y = stageStepLin.y.Data(:,1);
% u = stageStepLin.u.Data(:,1);

nd = 11;
Ts = 40e-6;
y = y(nd+1:end);
u = u(1:end-nd);

size(y)
size(u)


[A,B,C,D,K,R,~,s] = subid(y, u, 190, 12);

sys = ss(A, B, C, D*0, Ts);
sys.InputDelay = nd-1;

z = zero(sys);
gg = zpk([], z(1:2), 1, Ts);
sys = minreal(ss(gg)*sys/dcgain(gg), .01);
[A, B, C] = ssdata(sys);

F1 = figure(1); clf
frfBode(Gx_frf, w_s_x/2/pi, F1, '-k', 'Hz')


[~,~,ws] = bode(sys, {w_s_x(1), w_s_x(end)});

frfBode(sys, ws/2/pi, F1, '--r', 'Hz')
plotPZ_freqs(sys, F1, 0)
isstable(sys)
%
% figure(2)
% plot(s)
% 
% figure(3); clf
% pzplot(sys)
% 
% figure(4)
% step(sys)

% pzplot(sys)
% xlim([-2,2])
%
Qw = K*R*K';
L = dlqr(A', C', Qw, R)';

syscl = ss(A-L*C, B, C, 0, 40e-6);
figure(3); clf
pzplot(syscl, 'r')
hold on
pzplot(sys, 'b')
legend('closed loop', 'open loop')
%%
Fs = 1/Ts;
Ws = Fs*2*pi;

U = fft(u);
Y = fft(y);

N = length(U);
U = U(1:N/2);
Y = Y(1:N/2);
figure(1);
subplot(2,1,1)
freqs = [0:N/2-1]'/(N/2-1);
freqs = freqs*Fs/2;
% *Fs/(N);
% *(pi/(N/2));
% freqs = freqs*Ws/pi;
% freqs = freqs/2/pi;
semilogx(freqs, 20*log10(abs(Y./U)))

%%

modelFit.sys12 = sys;
modelFit.Gamma_w = K;
modelFit.Qw = Qw;


save(modelFitPath, 'modelFit')



