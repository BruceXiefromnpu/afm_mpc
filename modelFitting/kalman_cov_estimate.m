
clear;
clc

PHI = [0.75, -1.74, -0.3, 0, -0.15;
      0.09,   0.91, -0.0015, 0, -0.008;
      0,      0,     0.95,   0,  0;
      0,      0,      0,     0.55, 0;
      0,      0,      0,     0.0,  0.905;];
Gam = [0, 0, 0;  
       0, 0, 0;  
       24.64, 0, 0;
       0, 0.835, 0;
       0, 0, 1.83];
H = [1, 0, 0, 0, 1;
    0, 1, 0, 1, 0];

Q = eye(3);
R = eye(2);

Qo = diag([0.25, 0.5, 0.75]);
Ro = diag([0.4, 0.6]);

% True gain, M
M = dare(PHI', H', Gam*Q*Gam', R);
K = M*H'*inv(H*M*H' + R);
abs(eig(PHI - PHI*K*H))
%%
% Initial S.S. kalman gain:
Mo = dare(PHI', H', Gam*Qo*Gam', Ro);
Ko = Mo*H'*inv(H*Mo*H' + Ro);

N = 950;

ns = 5;
nu = 3;
no = 2;

x0 = zeros(ns,1);
xk = x0;
xk_hat = xk;

nu_vec = [];
for k=0:N
uk = mvnrnd(zeros(3,1), Q)';
vk = mvnrnd(zeros(2,1), R)';

% Actual dynamics simulation
xk_plus1 = PHI*xk + Gam*uk;

zk = H*xk + vk;

xk = xk_plus1;


% Estimate:
zk_hat = H*xk_hat;
xk_hat = xk_hat + Ko*(zk - zk_hat);
% Innovations
nu_k = zk - zk_hat;

% Propagate:
xk_hat = PHI*xk_hat;


% store:
nu_vec(1:2,k+1) = nu_k;
end

% -------------------------
% Compute estimates of Ck
Ck_mat=[]
Chat_stack = []
for k=0:N
    Ck = zeros(2,2);
    for i=k:N
        Ck = Ck + nu_vec(:,i+1)*nu_vec(:,i-k+1)';
    end
    Ck = Ck/size(nu_vec,2);
    Chat_stack = [Chat_stack; Ck];
    Ck_mat(:,:,k+1) = Ck;
end

BB = obsv(PHI, H)*PHI;
BB_sharp = (BB'*BB)\BB';

A=[];
for k=0:ns-1
    A = [A;
        H*(PHI*(eye(ns)-Ko*H))^k];
end
A = A*PHI;
%%
A_sharp = (A'*A)\A';
%%
clc
Co = inv(Ck_mat(:,:,1));
Chat_1_ns = Chat_stack(no+1:ns*no+no, :);

K_hat = Ko + A_sharp*Chat_1_ns*Co;

MH_hat = Ko*Co + A_sharp*Chat_1_ns;
H*MH_hat

R_hat = Co - H*MH_hat
%%

for k = 1:10
aa = PHI*(eye(ns) - K_hat*H);
qq = -PHI*(K_hat-Ko)*Co*(K_hat-Ko)'*PHI';
deltaM = dlyap(aa, qq);


MH_hat = MH_hat + deltaM*H';
Ko = K_hat;
K_hat = MH_hat*inv(H*MH_hat + R_hat);

norm(K_hat-Ko)
norm(K_hat)

end

abs(eig(PHI - PHI*K_hat*H))


