
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
B = Gam(:,1:2);
Q = eye(3);
R = eye(2);

QQ = Gam*Q*Gam';


% True gain, M
M = dare(PHI', H', QQ, R);
K = M*H'*inv(H*M*H' + R);
abs(eig(PHI - PHI*K*H))
%

N = 3950;

ns = 5;
nu = 3;
no = 2;

x0 = zeros(ns,1);
xk = x0;
xk_hat = xk;

nu_vec = [];
for k=0:N
uk = mvnrnd(zeros(5,1), QQ)';
vk = mvnrnd(zeros(2,1), R)';

% Actual dynamics simulation
xk_plus1 = PHI*xk + uk + B*[sin(k/N); cos(k/N)];

zk = H*xk + vk;

xk = xk_plus1;

% Propagate:
xk_hat = PHI*xk_hat + B*[sin(k/N); cos(k/N)];
y1 = H*xk_hat;

% store:
y2_vec(:,k+1) = zk - y1;

end

%%
% -------------------------
% Compute estimates of Ck
clc
p = 1000;
NN = 2900;
yy = y2_vec(:,1:NN);

% YY = zeros(no*(NN-p), NN);
YY = [];
for k = 1:p
    YY = [YY; 
        [zeros(no,k-1), yy(:,1:NN-k)]];
end
yp = yy(:,2:end);
alp_s = yy(:,2:end)/YY;

pp = alp_s*YY;
% alp_s = [eye(2), alp_s];

alps = reshape(alp_s, [], 2, p);


V = []
for k=0:p-1
   V = [V;H*PHI^k]; 
end


beta_s = [];
Y_K = [];
for k=1:p
    
   beta_k=alps(:,:,k);
   for i=1:k-1
      beta_k = beta_k + beta_s(:,:, k-i)*alps(:,:,i); 
   end
    beta_s(:,:,k) = beta_k;
   Y_K = [Y_K; beta_k];
   
end


% K_est = inv(V'*V)*V'*Y_K
K_est = V\Y_K
K

clc

abs(eig(PHI - K_est*H))
abs(eig(PHI - K*H))









