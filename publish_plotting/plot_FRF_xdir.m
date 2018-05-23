% The script will automatically save the resulting FRF's into a .mat file
% with the same name as sinesOut_FileName, replacing only .csv with .mat.
clear
clc

% plotting options
lw = 1.3;
saveOn = true; % set this to 0 to avoid saving the computed FRFs.


rmpath('C:\Users\arnold\Documents\MATLAB\miscScripts\system_id\')
addpath('functions')

addpath('/home/arnold/gradschool/sysID/matlab/functions')



% FC_data_file = 'x-axis_sines_infoFourierCoef_5-14-2018-03.csv';
FC_data_file = 'x-axis_sines_info_HIRESFourierCoef_5-14-2018-02.csv';

dataRoot = '/media/labserver/mpc-journal/sysID/FRF_data/'
FC_path = fullfile(dataRoot, FC_data_file);
load(fullfile(PATHS.sysid, 'FRF_data_current_stage2.mat'));

% For saving results:
frf_File = strrep(FC_path, '.csv', '.mat');

ssOpts = sweptSinesMeta('read', FC_path);

Ts = ssOpts.Ts;

[FC_s_raw, E_s, freqs, ssOpts] = SweptSines.read_FC_data(FC_path, ssOpts);

idx_uz = 1;
idx_stage = 2;
idx_powV = 3;
idx_powI = 4;

% Gains for current sensing 
R_sense = 0.1;
R3 = 68.4e3;
R4 = 10.15e6;
I_gain = (1/R_sense)*(R3/(R3+R4));
I_gain = 1/15.15; % Measured
% I_gain = 2.5757;
% Gains for amplifier output voltage

R2 = 1.732e6;
R1 = 29.7e6;
Vdiv_gain = R2/(R1+R2);
%Scale the Fourier coefficients of powV  by 1/Vdiv_gain and of powI by
%1/Rsense.
FC_s = FC_s_raw;
FC_s(:,idx_powI) = FC_s_raw(:, idx_powI)*I_gain;
FC_s(:,idx_powV) = FC_s_raw(:, idx_powV)/Vdiv_gain;

% Transfer functions
G_uz2stage = FC_s(:, idx_stage)./FC_s(:, idx_uz);
G_uz2powV = FC_s(:,idx_powV)./FC_s(:,idx_uz);
G_uz2powI = FC_s(:,idx_powI)./FC_s(:,idx_uz);
G_powV2powI = FC_s(:,idx_powI)./FC_s(:, idx_powV);

% Residual Energies

E_res_stage = (E_s(:, idx_stage) - 2*abs(0.5*FC_s_raw(:, idx_stage)).^2 )./E_s(:, idx_stage);
E_res_powI = (E_s(:, idx_powI) - 2*abs(0.5*FC_s_raw(:, idx_powI)).^2 )./E_s(:, idx_powI);
E_res_powV = (E_s(:, idx_powV) - 2*abs(0.5*FC_s_raw(:, idx_powV)).^2 )./E_s(:, idx_powV);

F2 = figure(2); clf
F2.PaperPosition = [1.3376    2.3454    5.8247    6.3093];

h1 = frfBode(G_uz2stage, freqs, F2, '-r', 'Hz');
h1.DisplayName = '$G_{y_x,u_x}$';
subplot(2,1,1)
set(gca, 'Position', [0.1300 0.5357 0.7750 0.3893]);
ylim([-100, 15])
title('Control signal to stage output')

yyaxis 'right'
h2 = semilogx(freqs, E_res_stage, ':k', 'LineWidth', 1.5);
h2.DisplayName = 'Normalized Residual Energy';
xlim([freqs(1), freqs(end)])
ylim([-0.1, 1.1])
leg1 = legend([h1, h2]);
set(leg1, 'Position', [0.1552    0.6722    0.3396    0.0637]);
subplot(2,1,2)
set(gca, 'Position', [0.1300 0.1100 0.7750 0.3543])
% saveas(F2, fullfile(PATHS.jfig, 'G_uz2stage_Eres.svg'))


% ----------------------------------------------------------------------- %
% Find the TF-bound

nm1 = modelFit.models.g_deluz2pow_1norm;
derz = tf([1 -1], 1, Ts)*nm1;
derz_frf = squeeze(freqresp(derz, freqs*2*pi));

% derz = tf([1 -1], 1, Ts);
% derz_frf = squeeze(freqresp(derz, freqs*2*pi));
% Derz_Guz2I = G_uz2powI./derz_frf;
% mag_max = max(abs(Derz_Guz2I));
% derz = mag_max*derz;
% derz_frf = derz_frf*mag_max;
% I_max = 0.1; %Amps
% fprintf('Mag-max = %.3f, deltaUk_max = %.3f\n', mag_max, I_max/mag_max);
     
      
F3 = figure(3); clf
ax3 = gca();
% F3.PaperPosition = [1.3376    2.3454    5.8247    6.3093];
h3 = semilogx(ax3, freqs, 20*log10(abs(G_uz2powV)), '-r', 'LineWidth', lw);
hold on, grid on;
h3.DisplayName = '$G_{V_X, u_X}$';

h4 = semilogx(ax3, freqs, 20*log10(abs(G_uz2powI)), '-k', 'LineWidth', lw);
h4.DisplayName = '$G_{I_x,u_x}$';

h5 = semilogx(ax3, freqs, 20*log10(abs(G_powV2powI)), '-b', 'LineWidth', lw);
h5.DisplayName = '$G_{I_x, V_X}$';

h10_bnd = semilogx(freqs, 20*log10(abs(derz_frf)),...
          ':k', 'LineWidth', lw);
h10_bnd.DisplayName = '$(z-1) ||g_{I_X,\Delta u_X}||_1$ ';
xlim([freqs(1), freqs(end)])

leg = legend([h3,h4,h5, h10_bnd]);
set(leg, 'location', 'SouthEast', 'FontSize', 14)
xlabel('Frequency [Hz]')
ylabel('Mag [dB]')

f4 = figure(4); clf
semilogx(freqs, E_res_powI);
hold on
semilogx(freqs, E_res_powV)
xlim([0.1, 10^4])
grid on
title('Residual energy')



if saveOn
  saveas(F2, fullfile(PATHS.jfig, 'G_uz2stage_Eres.svg'));
  saveas(F3, fullfile(PATHS.jfig, 'G_pow_and_current.svg'));
end
%%

modelFit.frf.G_uz2powV = G_uz2powV;
modelFit.frf.G_uz2stage = G_uz2stage;
modelFit.frf.G_uz2powI  = G_uz2powI;
modelFit.frf.G_powV2powI = G_powV2powI;
modelFit.frf.w_s       = freqs*2*pi;
modelFit.E_s = E_s;
modelFit.E_res_stage = E_res_stage;

modelFit.frf.freqs_Hz  = freqs;
modelFit.frf.Ts        = Ts;
modelFit.frf.freq_s    = freqs;



if saveOn
        save(frf_File, 'modelFit');
end










