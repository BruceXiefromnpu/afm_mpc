clear
% clc

saveon = true;
hash_opts = struct('Method', 'MD5', 'Format', 'hex', 'Input', 'file');

% load('hystID_data_4-30-2018_01.mat')
addpath(fullfile(PATHS.step_exp(), 'functions'))
addpath(fullfile(PATHS.step_exp(), 'functions', 'canon'))

[plants, ~, MF] = CanonPlants.plants_ns14(9,2);


% hyst_file = 'hystID_data_27-Aug-2018_01.mat';
% hyst_path = fullfile(PATHS.sysid, 'hysteresis', hyst_file);

hyst_path = MF.heritage.fpath_hystID;
% make sure the data hash not changed
hyst_hash_expected = MF.heritage.hyst_data_hash;
hyst_hash_current = DataHash(hyst_path, MF.heritage.hash_opts);
assert(strcmp(hyst_hash_expected, hyst_hash_current))

load(hyst_path)

% % To remove the vibrational and drift aspects, we need those models. Load them:
% modelFit_file = fullfile(PATHS.sysid, 'FRF_data_current_stage2.mat');
% load(modelFit_file)


Gvib = plants.Gvib;
Ts = Gvib.Ts;
gdrift = plants.gdrift;

ms = 1e3;

ux = hystData.u_exp;
yx = hystData.y_exp - mean(hystData.y_exp(1:500));
tvec = hystData.t_exp;

% Simulate the fitted hysteresis response fit
H = plants.hyst_sat;
u_hystsat = PIHyst.hyst_play_sat_op(ux, H.r, H.w, H.d, H.ws, H.r*0);
y_fit = lsim(Gvib*gdrift, u_hystsat, tvec);


% ux = downsample(ux, 50);
width = 3.5;
height = 3;
f1 = mkfig(1, width, height); clf

wd = 3.5;
ht = 3;

ax1 = axes('Units', 'inches', 'Position', [0.4550 0.35 3 2.58]);
h1 = plot(tvec, ux*dcgain(Gvib*gdrift), ':k');
h1.DisplayName = '$u_k G_d(0) G_{\textrm{vib}}(0)$';
hold on
h2 = plot(tvec, yx, 'k');
h2.DisplayName = 'Stage Response';

h3 = plot(tvec, y_fit, '--r');
h3.DisplayName = 'Model';
ylim([-8, 8]);
grid on
xlim([tvec(1), tvec(end)]);
leg1 = legend([h1, h2, h3]);
set(leg1, 'Position', [0.5546 0.7991 0.4325 0.1674], 'box', 'off');
xlab = xlabel('time [s]');
ylab = ylabel('Amplitude [$\mu$m]');

set(xlab, 'Units', 'inches', 'Position', [1.5, -0.18, 0])

el1 = annotation('ellipse', [0.6115 0.2696 0.0400 0.0500], 'Color', 'k');

ax2 = axes('Units', 'inches', 'Position', [2.5959 0.3583 0.86000 0.7000]);

a1 = annotation('line', [0.6256 0.7232], [0.3232 0.3715]); %, 'Color', 'k');
a2 = annotation('line', [0.631 0.7232], [0.2674 0.1389]); %, 'Color', 'k');
f1.CurrentAxes = ax2;
hold on, grid on
idx_start = 232080-1000;
idx_end = 240012+500;

plot(tvec(idx_start:idx_end), ux(idx_start:idx_end)*dcgain(Gvib*gdrift), ':k')
plot(tvec(idx_start:idx_end), yx(idx_start:idx_end), 'k')
plot(tvec(idx_start:idx_end), y_fit(idx_start:idx_end), '--r')

ylim([-5, -4])
xlim(Ts*[idx_start, idx_end])
set(ax2, 'XTick', [], 'YTick', [-4.5, -4], 'Box', 'on')




if saveon
  saveas(f1, fullfile(PATHS.jfig, 'hyst_response_ol.svg'))
end