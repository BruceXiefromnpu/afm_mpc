%%
rng(5);
saveon = true;
Ts = StageParams.Ts;
n_space = 2500;
N_repeat = 1;
N_steps = 20;
ymax = 7;
imps = -ymax + (2*ymax)*rand(N_steps,1);
imps = [imps;0; 7.0; -7.0; 0];

F1 = figure(1000); clf
step_ref = StepRef(imps, n_space);

step_ref.plot(F1);

step_ref.plot_settle_boundary(F1, 0.01, 'rel')


if saveon
  root = fullfile(PATHS.exp, 'step-exps');
% root = pwd;
  save(fullfile(root, 'many_steps_data_rand_ymax7.mat'), 'step_ref' );
end



