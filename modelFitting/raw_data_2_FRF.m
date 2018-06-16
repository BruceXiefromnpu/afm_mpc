% After doing the swept sines experiement in LabView, you should run this
% script. It's overall function is to process the sines of the reference
% input, control output, and system output and generate Fourier
% coefficients (by numerically integrating each signal against the base Fourier
% kernal for that frequency). The ratio of these coeffiecients will be the
% different frequency responses (roughly. In actuallity, we do some sneaky
% stuff with all the averages to help with noise). 

% The script will automatically save the resulting FRF's into a .mat file
% with the same name as sinesOut_FileName, replacing only .csv with .mat.

if ispc
  rmpath('C:\Users\arnold\Documents\MATLAB\miscScripts\system_id\')
  addpath('C:\Users\arnold\Documents\labview\sysID\matlab\functions')
else
  rmpath(fullfile(getMatPath(), 'toolboxes', 'system_id'))
  addpath('/home/arnold/gradschool/sysID/matlab/functions')
end

% addpath('functions')

clear
clc


% Experiemental data out file. Although we need the data from the
% input-data csv file, labview saves that file into the header of the
% output file. Labview will also save the data from the z-axis-settings
% cluster into the header. The following variable is the filename displayed
% in 'output-data-file-name' indicator in play_sysID_Z_Axis.vi'.


% FC_data_file = 'x-axis_sines_infoFourierCoef_4-16-2018stage-01.csv';
FC_data_file = 'x-axis_sines_infoFourierCoef_5-30-2018-01.csv';
dataRoot = fullfile(PATHS.sysid, 'FRF_data'); 
FC_path = fullfile(dataRoot, FC_data_file);

% For saving results:
frf_FileName = strrep(FC_data_file, '.csv', '.mat');
% save FRF data to .mat file?
save_data = true; 
% save the frf figure?
save_fig = false;


ssOpts = sweptSinesMeta('read', FC_path);

Ts = ssOpts.Ts;

[FC_s, E_s, freqs, ssOpts] = SweptSines.read_FC_data(FC_path, ssOpts);

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
FC_s(:,idx_powI) = FC_s(:, idx_powI)*I_gain;
FC_s(:,idx_powV) = FC_s(:, idx_powV)/Vdiv_gain;

G_uz2stage = FC_s(:, idx_stage)./FC_s(:, idx_uz);
G_uz2powV = FC_s(:,idx_powV)./FC_s(:,idx_uz);
G_uz2powI = FC_s(:,idx_powI)./FC_s(:,idx_uz);
G_powV2powI = FC_s(:,idx_powI)./FC_s(:, idx_powV);

F2 = figure(2); clf
F2.PaperPosition = [1.3376    2.3454    5.8247    6.3093];
h1 = frfBode(G_uz2stage, freqs, F2, '-r', 'Hz');
subplot(2,1,1)
title('Control signal to stage output')


F3 = figure(3); clf
F3.PaperPosition = [1.3376    2.3454    5.8247    6.3093];
h1 = frfBode(G_uz2powV, freqs, F3, '-r', 'Hz');
subplot(2,1,1)
title('Control signal to Amplifier Voltage');

F5 = figure(5); %clf
frfBode(G_uz2powI, freqs, F5, '--k', 'Hz');
frfBode(G_powV2powI, freqs, F5, '--m', 'Hz');

subplot(2,1,1)
title('Amplifier Current')
legend('From control signal', 'From powV');

f4 = figure(4);
semilogx(freqs, E_s(:, idx_powV));
hold on
semilogx(freqs, E_s(:, idx_powI));
xlim([0.1, 10^4])
grid on
title('Residual energy')




if save_fig
  saveas(F5, '../figures/pow-amp/G_uz2Current.svg');
end


modelFit.frf.G_uz2powV = G_uz2powV;
modelFit.frf.G_uz2stage = G_uz2stage;
modelFit.frf.G_uz2powI  = G_uz2powI;
modelFit.frf.G_powV2powI = G_powV2powI;
modelFit.frf.w_s       = freqs*2*pi;
modelFit.frf.freqs_Hz  = freqs;
modelFit.frf.Ts        = Ts;
modelFit.frf.freq_s    = freqs;



if save_data
        save(strrep(FC_path, '.csv', '.mat'), 'modelFit');
end










