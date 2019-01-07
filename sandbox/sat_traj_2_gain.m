clc
syscl = ss(sys_sim.a - sys_sim.b*K_lqr, sys_sim.b, sys_sim.c, 0, Ts);
% [~,~,X] = initial(syscl, Nx_sim);

N = 50;
X = zeros(4, N);
X(:,1) = -6*Nx_sim;
U = zeros(1, N);
for k=1:N
    
    uk = -K_lqr*X(:,k);
    
    uk = max(uk, -0.5);
    uk = min(uk, 0.5);
    
    X(:, k+1) = sys_sim.a*X(:,k) + sys_sim.b*uk;
    U(1, k) = uk;
end
X = X(:,1:end-1);
Y = sys_sim.C*X;
figure(1000)
subplot(2,1,1)
plot(Y)
subplot(2,1,2)
plot(U)

sat_idx = find(U < 0.5 & U>-0.5)
sat_idx(end) = [];
clc
Delta = U + K_lqr*X;

X0 = X(:, 1:end-1);
X1 = X(:, 2:end);
X1_min_AX0 = (X1 - sys_sim.a*X0)';


I = eye(length(Delta)-1);
I(sat_idx, sat_idx) = 0;

LHS = [I, -X0'];



DKT = LHS\(X1_min_AX0*pinv(sys_sim.b'));

K = DKT(end-3:end)'
K_lqr

% Do with a quadratic program. In the under-determined case above with
% ldiv, we get the minimum norm solution. SO the only difference here is
% that I make the "norm" not weight the K_lqr section.
H = blkdiag(I, eye(4)*0);
DKT = quadprog(H, [], [],[],LHS, (X1_min_AX0*pinv(sys_sim.b')));

KK = DKT(end-3:end)'
K_lqr


%%
U = step_data_clqr.results.opt_trajs_cell{1}.U_vec_s{25};
X = step_data_clqr.results.opt_trajs_cell{1}.X_vec_s{25};

U = U.Data(1:100)';
X = X.Data(1:100, :)';

sat_idx = find(U < 0.5 & U>-0.5)
sat_idx(end) = [];

X0 = X(:, 1:end-1);
X1 = X(:, 2:end);
X1_min_AX0 = (X1 - sys_sim.a*X0)';


I = eye(length(Delta)-1);
I(sat_idx, sat_idx) = 0;

LHS = [I, -X0'];

DKT = LHS\(X1_min_AX0*pinv(sys_sim.b'));

K = DKT(end-3:end)'
K_lqr


X = zeros(4, N);
X(:,1) = -6.1*Nx_sim;
U = zeros(1, N);

for k=1:N
    
    uk = -K_lqr*X(:,k);
    
    uk = max(uk, -0.5);
    uk = min(uk, 0.5);
    
    X(:, k+1) = sys_sim.a*X(:,k) + sys_sim.b*uk;
    U(1, k) = uk;

