%%

n_space = 800;
N_repeat = 3;
imps = [ 1, 2.5,-2,1.25, -2, 1, -3, -.5]';
imps = [imps; -sum(imps)];
imps = [0; repmat(imps, N_repeat, 1)];

ref_s = cumsum(imps)
N_imp = length(imps)


impulse_idx= (1:n_space:N_imp*n_space)';
u_vec = zeros((N_imp+1)*(n_space), 1);
u_vec(impulse_idx) = imps;
u_vec = cumsum(u_vec);
% u_vec = repmat(cumsum(u_vec), 3,1);

t_vec = (0:length(u_vec)-1)'*Ts;
figure
plot(t_vec, u_vec);
grid on
%%
clc
ref_traj_params.ref_traj = timeseries(u_vec(:), t_vec(:));
ref_traj_params.impulse_idx = impulse_idx;
ref_traj_params.ref_s = ref_s;
save('many_steps.mat', 'ref_traj_params' )
