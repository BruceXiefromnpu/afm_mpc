%%
Ts = StageParams.Ts;
n_space = 1600;
N_repeat = 3;
imps = [ 1, 2.5,-2,1.25, -2, 1, -3, -.5]';
imps = [imps; -sum(imps)];
imps = [0; repmat(imps, N_repeat, 1)];

step_ref = StepRef(imps, n_space);


root = fullfile(PATHS.exp, 'step-exps', 'many_steps_data');
save(fullfile(root, 'many_steps.mat'), 'step_ref' );


