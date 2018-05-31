%%
rng(5);
saveon = true;
Ts = StageParams.Ts;
n_space = 2500;
N_repeat = 1;
N_steps = 20;
ymax = 7;
imps = -ymax + (2*ymax)*rand(N_steps,1)
% imps = [ 1, 2.5,-2,1.25, -2, 1, -3, -.5]';
imps = diff([0; imps]);
imps = [imps; -sum(imps); 7.0; -14; 7.0];
imps = [0; repmat(imps, N_repeat, 1)];

ref_s = cumsum(imps)
N_imp = length(imps)


impulse_idx= (1:n_space:N_imp*n_space)';
u_vec = zeros((N_imp+1)*(n_space), 1);
u_vec(impulse_idx) = imps;
u_vec = cumsum(u_vec);
% u_vec = repmat(cumsum(u_vec), 3,1);

t_vec = (0:length(u_vec)-1)'*Ts;
figure(1); clf
plot(t_vec, u_vec);
grid on


ref_traj_params.ref_traj = timeseries(u_vec(:), t_vec(:));
ref_traj_params.impulse_idx = impulse_idx;
ref_traj_params.ref_s = ref_s;
if saveon
%   root = fullfile(PATHS.exp, 'step-exps', 'many_steps_data_rand');
root = pwd;
  save(fullfile(root, 'many_steps_data_rand_ymax7.mat'), 'ref_traj_params' );
end



