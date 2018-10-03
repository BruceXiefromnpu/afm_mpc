clc
exp_path = fullfile(PATHS.step_exp, 'many_steps_sweep_gamma_21-Sep-2018_01');
ts_total_file = 'ts_totalconst-sig_09-21-2018.mat';
fpath = fullfile(exp_path, ts_total_file);
dat_CS = load(fpath)

dat_CS.ts_sum_sim_results(:,2:3) = dat_CS.ts_sum_sim_results(:,2:3)*1000;
dat_CS.ts_sum_exp_results(:,2:3,:) = dat_CS.ts_sum_exp_results(:,2:3,:)*1000;

ts_sum_exp_results_CS = dat_CS.ts_sum_exp_results(:,2:end,:);
ts_exp_means_CS = mean(ts_sum_exp_results_CS, 3); % average across pages (3rd dim)
ts_exp_stdv_CS = sqrt(var(ts_sum_exp_results_CS,0,3)); % average across pages (3rd dim)
% ts_data_CS = struct('means', ts_exp_means_CS, 'stdv', ts_exp_stdv_CS, 'gams', gams_ts);


S = '$\gamma$&MPC-sim & SLF-sim & MPC-exp & SLF-exp\\';
S = sprintf('%s\n\\toprule\n', S)

for k=1:length(dat_CS.gam_s)
  gam_k = dat_CZ.gam_s(k);
  if gam_k <= 0.1
    gam_str = sprintf('$10^{%d}$', log10(gam_k));
  else
    gam_str = sprintf('%g', gam_k);
  end
  
  if isnan(dat_CS.ts_sum_sim_results(k,2))
  S = sprintf('%s%s&%.1f&--&%.1f&--\\\\\n', S,gam_str,...
    dat_CS.ts_sum_sim_results(k,3), ts_exp_means_CS(k,2));
  else
  S = sprintf('%s%s&%.1f&%.1f&%.1f&%.1f\\\\\n', S, gam_str,...
    dat_CS.ts_sum_sim_results(k,3), dat_CS.ts_sum_sim_results(k,2),...
    ts_exp_means_CS(k,2), ts_exp_means_CS(k,1));
  end
  
  
end
S
CS_tspath = fullfile(PATHS.MPCJ_root, 'latex', 'ts_means_CS.tex')
fid = fopen(CS_tspath, 'w+');
fprintf(fid, '%s', S);
fclose(fid);




%%
exp_path = fullfile(PATHS.step_exp, 'many_steps_sweep_gamma_21-Sep-2018_01');
ts_total_file = 'ts_totalchoose-zet_09-21-2018.mat';
fpath = fullfile(exp_path, ts_total_file);
dat_CZ = load(fpath)

dat_CZ.ts_sum_sim_results(:,2:3) = dat_CZ.ts_sum_sim_results(:,2:3)*1000;
dat_CZ.ts_sum_exp_results(:,2:3,:) = dat_CZ.ts_sum_exp_results(:,2:3,:)*1000;

ts_sum_exp_results_CZ = dat_CZ.ts_sum_exp_results(:,2:end,:);
ts_exp_means_CZ = mean(ts_sum_exp_results_CZ, 3); % average across pages (3rd dim)
ts_exp_stdv_CZ = sqrt(var(ts_sum_exp_results_CZ,0,3)); % average across pages (3rd dim)
% ts_data_CS = struct('means', ts_exp_means_CS, 'stdv', ts_exp_stdv_CS, 'gams', gams_ts);


S = '$\gamma$&MPC-sim & SLF-sim & MPC-exp & SLF-exp\\';
S = sprintf('%s\n\\toprule\n', S)

for k=1:length(dat_CZ.gam_s)
  gam_k = dat_CZ.gam_s(k);
  if gam_k <= 0.1
    gam_str = sprintf('$10^{%d}$', log10(gam_k));
  else
    gam_str = sprintf('%g', gam_k);
  end
  if isnan(dat_CZ.ts_sum_sim_results(k,2))
    S = sprintf('%s%s&%.1f&--&%.1f&--\\\\\n', S, gam_str,...
      dat_CZ.ts_sum_sim_results(k,3), ts_exp_means_CZ(k,2));
  else
    S = sprintf('%s%s&%.1f&%.1f&%.1f&%.1f\\\\\n', S, gam_str,...
      dat_CZ.ts_sum_sim_results(k,3), dat_CZ.ts_sum_sim_results(k,2),...
      ts_exp_means_CZ(k,2), ts_exp_means_CZ(k,1));
  end
  
end
S
CS_tspath = fullfile(PATHS.MPCJ_root, 'latex', 'ts_means_CZ.tex')
fid = fopen(CS_tspath, 'w+');
fprintf(fid, '%s', S);
fclose(fid);


% function str = gam_str()
%   
%   
%   
% end
