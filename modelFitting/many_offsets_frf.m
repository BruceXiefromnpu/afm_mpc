clc
addpath('/home/arnold/gradschool/sysID/matlab/functions/');
data_root = '/home/arnold/matlab/afm_mpc_journal/data/sys_id/many_offsets';

data_names = ls(data_root);
data_names = strsplit(strip(data_names), '\n')'


four_idx = ~cellfun(@isempty, strfind(data_names, 'Fourier'))

data_names_FC = data_names(four_idx')
data_names_time = data_names(~four_idx')


data_names_FC = data_names_FC(1)
data_names_time = data_names_time(1)
%%
FRF_cell = {};
%%
F1 = figure(1); clf
F2 = figure(2); clf
F3 = figure(3); clf
F4 = figure(4); clf
%%
h_uz2stage = gobjects(length(data_names_FC),1);
h_uz2powV = gobjects(length(data_names_FC),1);
h_uz2powI = gobjects(length(data_names_FC),1);
h_powV2powI = gobjects(length(data_names_FC),1);

LS = {'-r', '--b', ':k', '--m'};

for k = 1:length(data_names_FC)
  
  fpath = fullfile(data_root, data_names_FC{k});
  
  ssOpts = sweptSinesMeta('read', fpath);
  
  [FC_data, E_s, freqs, ssOpts] = SweptSines.read_FC_data(fpath, ssOpts);
  frf_dat_k.ssOpts = ssOpts;
  frf_dat_k.FC_data = FC_data;
  frf_dat_k.E_s = E_s;
  frf_dat_k.freqs = freqs;
  FRF_dat_cell{k} = frf_dat_k;
  
  G_uz2stage = FC_data(:,2)./FC_data(:,1);
  h_k = frfBode(G_uz2stage, freqs, F1, LS{k}, 'Hz');
  h_k.DisplayName = sprintf('$G_{uz2stage}$, offset=%.2f, Amp = %.1f', ssOpts.offset, ssOpts.Amps(1));
  h_uz2stage(k) = h_k;
  
  G_uz2powV = FC_data(:,3)./FC_data(:,1);
  h_k = frfBode(G_uz2powV, freqs, F2, LS{k}, 'Hz');
  h_k.DisplayName = sprintf('$G_{uz2powV}$, offset=%.2f, Amp = %.1f', ssOpts.offset, ssOpts.Amps(1));
  h_uz2powV(k) = h_k;
  
  G_uz2powI = FC_data(:,4)./FC_data(:,1);  
  h_k = frfBode(G_uz2powI, freqs, F3, LS{k}, 'Hz');
  h_k.DisplayName = sprintf('$G_{uz2powI}$, offset=%.2f, Amp = %.1f', ssOpts.offset, ssOpts.Amps(1));
  h_uz2powI(k) = h_k;
  
  G_powV2powI = FC_data(:,4)./FC_data(:,3);
  h_k = frfBode(G_powV2powI, freqs, F4, LS{k}, 'Hz');
  h_k.DisplayName = sprintf('$G_{powV2powI}$, offset=%.2f, Amp = %.1f', ssOpts.offset, ssOpts.Amps(1));
  h_powV2powI(k) = h_k;
  
  
end
%%
k= 1
  G_uz2stage = FC_data(:,2)./FC_data(:,1);
  h_k = frfBode(G_uz2stage, freqs, F1, LS{k}, 'Hz');
  h_k.DisplayName = sprintf('$G_{uz2stage}$, offset=%.2f, Amp = %.1f', ssOpts.offset, ssOpts.Amps(1));
  h_uz2stage(k) = h_k;
  
  G_uz2powV = FC_data(:,3)./FC_data(:,1);
  h_k = frfBode(G_uz2powV, freqs, F2, LS{k}, 'Hz');
  h_k.DisplayName = sprintf('$G_{uz2powV}$, offset=%.2f, Amp = %.1f', ssOpts.offset, ssOpts.Amps(1));
  h_uz2powV(k) = h_k;
  
  G_uz2powI = FC_data(:,4)./FC_data(:,1);  
  h_k = frfBode(G_uz2powI, freqs, F3, LS{k}, 'Hz');
  h_k.DisplayName = sprintf('$G_{uz2powI}$, offset=%.2f, Amp = %.1f', ssOpts.offset, ssOpts.Amps(1));
  h_uz2powI(k) = h_k;
  
  G_powV2powI = FC_data(:,4)./FC_data(:,3);
  h_k = frfBode(G_powV2powI, freqs, F4, LS{k}, 'Hz');
  h_k.DisplayName = sprintf('$G_{powV2powI}$, offset=%.2f, Amp = %.1f', ssOpts.offset, ssOpts.Amps(1));
  h_powV2powI(k) = h_k;
figure(F1)
legend(h_uz2stage);
subplot(2,1,1)
title('Uz 2 stage')

figure(F2)
legend(h_uz2powV);
subplot(2,1,1)
title('Uz 2 powV')


figure(F3)
legend(h_uz2powI);
subplot(2,1,1)
title('Uz 2 powI')

figure(F4)
legend(h_uz2powI);
subplot(2,1,1)
title('powV 2 powI')


%%

fname = fullfile(data_root, data_names_time{1})

ssOpts2 = sweptSinesMeta('read', fname)

chans_idx = [1, 2,3,4];
[xk, ssOpts2] = SweptSines.sys_id_raw_data(ssOpts2, chans_idx);
%%

ssOpts_tmp = ssOpts2;
% ssOpts_tmp.num_periodsMT = ssOpts_tmp.num_periodsMT.*ssOpts_tmp.NumAve;
% ssOpts_tmp.NumAve = ssOpts_tmp.NumAve*0+1;
xk_FC1_mat_page = rawData2FourierCoeffs(xk, ssOpts_tmp);

% Do averaging. Cross data is only accurate if we calculate it for each 
% measurement and then average the results (not average first and then 
% calc).

G_xdir = zeros(ssOpts2.num_channels, ssOpts2.num_channels, length(ssOpts2.freq_s));
for kk=1:ssOpts2.num_channels
  for jj = 1:ssOpts2.num_channels
    G_xdir(kk,jj, :) = avgFreqData(xk_FC1_mat_page(:,:,kk),...
                                       xk_FC1_mat_page(:,:,jj), ssOpts2);
  end
end

G_uz2stage = squeeze(G_xdir(2,1,1:end-1)./G_xdir(1,1,1:end-1));
G_uz2powV = squeeze(G_xdir(3,1,1:end-1)./G_xdir(1,1,1:end-1));
G_uz2powI = squeeze(G_xdir(4,1,1:end-1)./G_xdir(1,1,1:end-1));

frf_dat_k.ssOpts = ssOpts2;
frf_dat_k.CrosSpecData = G_xdir;
frf_dat_k.freqs = ssOpts2.freq_s_adjusted;
FRF_dat_cell_breakup{k} = frf_dat_k;
%%
F10 = figure(1); %clf
F12 = figure(2); %clf
F13 = figure(3); %clf
F14 = figure(4); %clf
k = 2;

h_k = frfBode(G_uz2stage, freqs, F10, LS{k}, 'Hz');
h_k.DisplayName = sprintf('$G_{uz2stage}$, offset=%.2f, Amp = %.1f', ssOpts2.offset, ssOpts.Amps(1));
h_uz2stage(k) = h_k;


h_k = frfBode(G_uz2powV, freqs, F12, LS{k}, 'Hz');
h_k.DisplayName = sprintf('$G_{uz2powV}$, offset=%.2f, Amp = %.1f', ssOpts2.offset, ssOpts.Amps(1));
h_uz2powV(k) = h_k;


h_k = frfBode(G_uz2powI, freqs, F13, LS{k}, 'Hz');
h_k.DisplayName = sprintf('$G_{uz2powI}$, offset=%.2f, Amp = %.1f', ssOpts2.offset, ssOpts.Amps(1));
h_uz2powI(k) = h_k;


h_k = frfBode(G_powV2powI, freqs, F14, LS{k}, 'Hz');
h_k.DisplayName = sprintf('$G_{powV2powI}$, offset=%.2f, Amp = %.1f', ssOpts2.offset, ssOpts.Amps(1));
h_powV2powI(k) = h_k;



figure(F10)
legend(h_uz2stage);
subplot(2,1,1)
title('Uz 2 stage')

figure(F12)
legend(h_uz2powV);
subplot(2,1,1)
title('Uz 2 powV')


figure(F13)
legend(h_uz2powI);
subplot(2,1,1)
title('Uz 2 powI')

figure(F14)
legend(h_uz2powI);
subplot(2,1,1)
title('powV 2 powI')

















