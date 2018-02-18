function [ A, B, C, D] = frf_id_nonunif_mckelvey(Gfrf, w_s,  R_s, opts)

q = opts.q;
n = opts.n;
p = opts.p;
doD = opts.Feedthrough;
Ts=1;
Gk_s = Gfrf;

M = length(w_s);
% k_s = w_s*M*Ts/pi;
% plot(k_s)
% Build the Big matrices. GG & Wm defined in eqs (47) & (48).
% Are q x M

[GG, Wm] = makeGG_WM(Gk_s, w_s, Ts, q, M);


% Noise:
alp = 1;
KKT = alp*real( Wm*diag(R_s)*Wm');
KKT = (KKT + KKT')/2;
% keyboard
% The paper defines K * K^T, but chol returns R s.t R^T * R = A. Thus, take
% the transpose.
% keyboard
% 
% KT = chol(KKT); 
% K = KT';
% cholesky seems to always be non-positive definite and thus fails. Since
% KKT should be symmetric, the svd should also give a square root, though
% won't be lower triangular.
[UK, SK, VK] = svd(KKT);
K = UK*sqrt(SK);
sk = diag(SK);
figure
plot(log(sk))

W_ = [real(Wm), imag(Wm)];
G_ = [real(GG), imag(GG)];

WG = [W_;
      G_];

% The paper says we should use QR factorization to get something equivalent
% to W_orthog. This allows us to avoid the inversion necessary to compute
% W_othog = eye(2*M) - W_'*inv(W_*W_')*W_. I had trouble with this, since
% we are supposed to end up with G_*W_orthog = R22*QT2, but R22 was always
% zero. Using the economing version fixes this.
%
% Also, we want the QR factorization such that
%
% [ W_ ] = [R11  0 ][Q1^T]
% [ G_ ]   [R21 R22][Q2^T]
% But qr(A) returns Q*R, with R upper triangular. Thus, we take the qr
% factorization of WG^T. 
% 
%  If A' = Q*R, then A = R'*Q', with R' lower triangular. 
%

[Q, RT] = qr(WG', 0);
R = RT';

% We want to partition Q^T between rows, so partition Q between columns.
Q1 = Q(:,1:q);
Q2 = Q(:,q+1:end);

% R11 = R(1:q, 1:q);
% R21 = R(q+1:end, 1:end-q);
R22 = R(q+1:end, end-q+1:end);

[U, S, V] = svd(K\R22);

% Plot singular values.
figure(3); clf;
plot(diag(S))
title('Singular values of K^{-1}*R_{22}')

[~, n_U] = size(U);

Us = U(:,1:n);

% Identity selectors
J1 = [eye( (q-1)*p) zeros((q-1)*p, p)];
J2 = [zeros((q-1)*p, p), eye((q-1)*p)];
J3 = [eye(p), zeros(p, (q-1)*p)];


X = J1*K*Us;
% Back out A & C
A = pinv(X)*J2*K*Us;
C = J3*K*Us;

% B & D via least squares:
% Least square problem to solve, for each w_k,
% G_k ~= [C^(exp(jw)I - A)^-1, 1][ B ]
%                                [ D ]
% For B and D. For each frequency, we stack these up. With covarince data
% included, this becomes:
% R^(-1/2)G_k ~= R^(-1/2)[C^(exp(jw)I - A)^-1, 1][ B ]
%                                                [ D ]
% For B and D. For each frequency, we stack these up

CIA = [];
GG = [];
for i = 1:M
    GG = [GG; R_s(i)^(-1/2)*Gk_s(i)];
%                             [C*inv(exp(1i*w_s(i)*Ts)*eye(n) - A)
    if doD
        CIA = [CIA; R_s(i)^(-1/2)*[C/(exp(1i*w_s(i)*Ts)*eye(n) - A), 1] ];
    else
        CIA = [CIA; R_s(i)^(-1/2)*[C/(exp(1i*w_s(i)*Ts)*eye(n) - A)] ];  
    end
    
end
% 
CIA = [real(CIA); imag(CIA)];
GG  = [real(GG); imag(GG)];
% keyboard
BD = lsqlin(CIA, GG);
B = BD(1:n);
if doD
    D = BD(n+1:end);
else
    D = 0;
end

end

function [GG, Wm] = makeGG_WM(Gk_s, wk_s, Ts, q, M)

    GG = zeros(q, M);
    Wm = zeros(q, M);
    for i = 1:M % across cols
        k_s = [0:1:q-1]';
        w_i = wk_s(i);
%         exp_jw_ks = exp(1i*w_i*Ts*k_s);
        exp_jw_ks = exp(1i*w_i*k_s); 
        GG(:,i) = exp_jw_ks*Gk_s(i);
        Wm(:,i) = exp_jw_ks;
    end

    GG = GG/sqrt(M);
    Wm = Wm/sqrt(M);
end













