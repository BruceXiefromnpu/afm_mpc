clear

addpath(fullfile(getMatPath, 'afm_mpc_journal', 'functions'));

exp_path = fullfile(PATHS.step_exp, 'many_steps_sweep_gamma10-Sep-2018_02');
ts_total_file = 'ts_totalconst-sig_09-10-2018.mat';

fpath = fullfile(exp_path, ts_total_file);
load(fpath)
whos

gams = ts_sum_exp_results(:,1,1)
ts_sum_exp_results = ts_sum_exp_results(:,2:end,:)


ts_exp_means = mean(ts_sum_exp_results, 3); % average across pages (3rd dim)
ts_exp_stdv = sqrt(var(ts_sum_exp_results,0,3)); % average across pages (3rd dim)

figure(1); clf

hlin_sim = semilogx(gams, ts_sum_sim_results(:,2)*1000, 'o');
hold on
hmpc_sim = semilogx(gams, ts_sum_sim_results(:,3)*1000, 'o');

% experimental

errorbar(gams, ts_exp_means(:,1)*1000, ts_exp_stdv(:,1)*1000, 'r.');
errorbar(gams, ts_exp_means(:,2)*1000, ts_exp_stdv(:,2)*1000, 'k.');

grid on


xlabel('$\gamma$')
ylabel('sum of settle times [ms]')