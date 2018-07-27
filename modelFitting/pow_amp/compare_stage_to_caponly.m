% Plot the stage vs simple cap circuit.
clear
clc

MF_stage = load('FRF_data_current_stage2.mat');
MF_cap   = load('FRF_data_current_caponly2.mat');

MF_stage = MF_stage.modelFit;
MF_cap = MF_cap.modelFit;

freqs_cap = MF_cap.frf.freq_s;
freqs_stage = MF_stage.frf.freq_s;

% Find the TF-bound
derz = tf([1 -1], 1, MF_stage.frf.Ts);
derz_frf = squeeze(freqresp(derz, freqs_stage*2*pi));
Derz_Guz2I = MF_stage.frf.G_uz2powI./derz_frf;
mag_max = max(abs(Derz_Guz2I));
derz = mag_max*derz;
derz_frf = derz_frf*mag_max;
I_max = 0.1; %Amps
fprintf('Mag-max = %.3f, deltaUk_max = %.3f\n', mag_max, I_max/ ...
        mag_max);


F10 = figure(10); clf
% F10.PaperPosition = [1.3333, 2.4635, 5.8333, 6.0729];
F10.Position = [-1833, 408, 560, 470];
lw = 1.5;
h10_1 = semilogx(freqs_stage, 20*log10(abs(MF_stage.frf.G_powV2powI)),...
                 '-r', 'LineWidth', lw);
h10_1.DisplayName = '$G_{I_X,V_X}$ (stage)';
hold on, grid on

h10_bnd = semilogx(freqs_stage, 20*log10(abs(derz_frf)),...
          ':k', 'LineWidth', lw);
h10_bnd.DisplayName = '$\tilde{G}_{I_X,u_X}$ ';

h10_2 = semilogx(freqs_cap, 20*log10(abs(MF_cap.frf.G_powV2powI)),...
  '--b', 'LineWidth', lw);
h10_2.DisplayName = '$G_{I_X,V_X}$ (dummy cap)';

h10_3 = semilogx(freqs_stage, 20*log10(abs(MF_stage.frf.G_uz2powI)),...
  '-k', 'LineWidth', lw);
h10_3.DisplayName = '$G_{I_X,u_X}$ (stage)';

h10_4 = semilogx(freqs_cap, 20*log10(abs(MF_cap.frf.G_uz2powI)),...
  '--g', 'LineWidth', lw);
h10_4.DisplayName = '$G_{I_X,u_X}$ (dummy cap)';

title('control Voltage to power-amp current');
leg1 = legend([h10_1, h10_2, h10_3, h10_4, h10_bnd]);
set(leg1, 'location', 'northwest', 'FontSize', 14)
xlim([freqs_stage(1), freqs_stage(end)]);
xlabel('frequency [Hz]')
ylabel('Mag [dB]')


F30 = figure(30); clf
F30.PaperPosition = [1.3333, 2.4635, 5.8333, 6.0729];
[h30_1m, h30_1p] = frfBode(MF_stage.frf.G_uz2stage, freqs_stage, F30, 'Hz', '-r');
[h30_2m, h30_2p] = frfBode(MF_stage.frf.G_uz2powV, freqs_stage, F30, 'Hz', '--b');
[h30_3m, h30_3p] = frfBode(MF_stage.frf.G_uz2powI, freqs_stage, F30, 'Hz', '-k');
[h30_4m, h30_4p] = frfBode(MF_stage.frf.G_powV2powI, freqs_stage, F30, 'Hz', ':k');

h30_1p.DisplayName = '$G_{y_X,u_X}$';
h30_2p.DisplayName = '$G_{V_X,u_X}$';
h30_3p.DisplayName = '$G_{I_X,u_X}$';
h30_4p.DisplayName = '$G_{I_X,V_X}$';

subplot(2,1,2)
leg1 = legend([h30_1p,h30_2p,h30_3p,h30_4p]);
set(leg1, 'location', 'southwest', 'FontSize', 14);
ylim([-1000, 300])

saveon = 0;

if saveon == 1
  saveas(F10, fullfile(PATHS.jfig, 'stage_vs_cap.svg'))
  saveas(F30, fullfile(PATHS.jfig, 'stage_frfs_all.svg'))
elseif saveon ==2
  saveas(F10, 'figures/stage_vs_cap.svg')
  saveas(F30, 'figures/stage_frfs_all.svg')
end
