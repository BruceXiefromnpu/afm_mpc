%%
Ts = StageParams.Ts;
n_space = 1600;
N_repeat = 3;
imps = [ 1, 2.5, -2.5, -2,1.25,.5, -2, 1, -3, -.5]';
imps = [imps; 0];
% imps = [0; repmat(imps, N_repeat, 1)];

step_ref = StepRef(imps, n_space);

step_ref.plot

root = fullfile(PATHS.exp, 'step-exps');
save(fullfile(root, 'many_steps.mat'), 'step_ref' );


