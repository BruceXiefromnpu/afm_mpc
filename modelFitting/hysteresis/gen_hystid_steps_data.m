
clear
Ts = 40e-6;
u_max = 7;
n_space = 2000;
n_up = 5;
step_sz = u_max/n_up;

% imps = [ones(9,1); -ones(18,1)];

imps = [0;step_sz*ones(n_up, 1); -step_sz*ones(2*n_up,1)];

for k = 2*n_up-1:-1:1
  if mod(k,2) == 0
    sgn = -1;
  else
    sgn = 1;
  end
  imps = [imps; sgn * step_sz * ones(k, 1)];
  
end

% ref_s = cumsum(imps)
N_imp = length(imps)


impulse_idx= (1:n_space:N_imp*n_space)';
u_vec = zeros((N_imp)*(n_space), 1);
u_vec(impulse_idx) = imps;
u_vec = cumsum(u_vec);
% u_vec = repmat(cumsum(u_vec), 3,1);

t_vec = (0:length(u_vec)-1)'*Ts;
figure(1); clf
plot(t_vec, u_vec);
grid on
%%
save('hystid_steps_input.m', 'u_vec', 'imps', 'u_max', 'n_space', 'step_sz',...
  'n_up', 'Ts');
