clc
[a, b] = ssdata(sys_sim);
R = 0;
ns = size(b,1);
nu = 1;
Q = eye(ns);
S = zeros(ns, 1);
N =23;
% -------------------------------------------------------------------------
% Sparse cLQR problem

% Build the cost function. THis should look like:
%      [R    |S 0 0 ][u0]
%      [   R |0 S 0 ][u1]
%  H = [-----|------][--]
%      [S' 0 |Q     ][x0]
%      [0  S'|  Q   ][x1]
%      [0  0 |    Qp][xN]

% First, build the part without the cross term. Use eye(N) because
% 0...N-1 has N elements.

RRR = kron(eye(N), R);
wghts = 0*linspace(1,1000, N+1);
QQQp = [];
for k=1:N+1
    QQQp = blkdiag(QQQp, eye(ns)*wghts(k));
end
% QQQp = blkdiag(kron(eye(N), Q), Qp);
Sx = [kron(eye(N), S);
  zeros(ns, N*nu)]; % The lower left.
Su = Sx';             % The upper right.

% The main event:
H = [RRR, Su;
    Sx, QQQp];


% Now build the equality constraint. This should look like:
%  [0 0  0|    |0  0 ][u0]   [x0]
%  [0 0  0| I  |0  0 ][u1]   [0]
%  [0 0  0|    |0  0 ][uN]   [0]
%  [-----------------][--] = [--]
%  [B 0  0 |A -I  0 ][x0]    [0 ]
%  [0 B  0 |0  A  -I][x1]    [0 ]
%                    [x_N] 

I_b = eye(N, N);
Aeq_b = kron(I_b, b);
Aeq_a = [kron(eye(N), a), zeros(N*ns, ns)]...
        + [zeros(N*ns, ns), kron(eye(N), -eye(ns))];


Aeq = [[zeros(ns, nu*N), eye(ns), zeros(ns, N*ns)];
       [Aeq_b, Aeq_a];
       zeros(ns, nu*N+ns+(N-1)*ns), eye(ns)];


beq = zeros(ns*(N+2), 1);

beq(1:ns) = -xf;
lb = repmat(-0.5, N,1);
ub = -lb;

UX = quadprog(H, [], [], [], Aeq, beq, lb, ub);

U = UX(1:N);
X = reshape(UX(N+1:end), 4, []);
Y = sys_sim.c*X;


U = [U; zeros(2*N, 1)];
[y, t, x] = lsim(sys_sim, U, [0:1:length(U)-1]'*Ts, -xf);

figure(6), plot(t, y)
figure(5), plot(t, U)




%%
% this is the convex relaxation idea. It doesn't seem to really work that
% well. I think I have this right...
% % % N = 50;
% % % P_tmp = path(); % so we can reset the path when we're done
% % % addpath(genpath(fullfile(getMatPath, 'solvers/cvx')))
% % % wghts = linspace(1,1000, N);
% % % wghts = log10(wghts);
% % % x00 = -xf;
% % % 
% % % cvx_begin
% % % variable u(N)
% % % % variable J
% % % J = 0;
% % % for k=1:N
% % %    delt=[];
% % % 
% % %     for j=0:k-1
% % %         delt = [a^j*b delt];
% % %     end
% % %     
% % %    J = J + log(k)*norm(a^k * x00 + delt*u(1:k), 2);
% % %     
% % % end
% % % 
% % % minimize J
% % % 
% % % subject to
% % % 
% % % norm(u, Inf) <= 0.5
% % % 
% % % 
% % % cvx_end
% % % 
% % % % path(P_tmp)
% % % %%
% % % 
% % % 
% % % uu = [u; zeros(2*N,1)];
% % % 
% % % [y, t, x] = lsim(sys_sim, uu, [0:1:length(uu)-1]'*Ts, x00);
% % % 
% % % figure(5)
% % % plot(t, uu)
% % % 
% % % figure(6)
% % % plot(t, y)


