clc

figure(1);clf
subplot(2,1,1)
plot(y)

subplot(2,1,2)
plot_ind(u)


kplus = [4,19]
% kmin  = 
xf = Nx_sim*ref_k;
xerr = -X.Data +repmat(xf', N, 1);

k_satp = find(u >= 0.5-.0001)
k_satm = find(u <= -0.5 +.0001)

k_usat = [5, 12, 20];

% x_usat = xf - X.Data(k_usat', :)
x_usat = xerr(k_usat', :)
u_usat  = U.Data(k_usat)

X_satp = xerr(k_satp, :);
X_satm = xerr(k_satm, :);
u_satp = u(k_satp);
u_satm = u(k_satm);

X_sat = [-X_satp;
         X_satm];
u_sat = [-repmat(0.5, length(k_satp), 1);
          repmat(-0.5, length(k_satm), 1)];
      

% K_est = (x_usat\u_usat)'
K_est = lsqlin(x_usat, u_usat, X_sat, u_sat)';


N = length(U.Data);

uvec_est = zeros(N, 1);
yvec_est = zeros(N, 1);
xvec_est = zeros(4, N);
xvec_est(:, 1) = xerr(1,:);


uvec_lqr = zeros(N, 1);
yvec_lqr = zeros(N, 1);
xvec_lqr = zeros(4, N);
xvec_lqr(:, 1) = xerr(1,:);


xk = xerr(1,:)'*0;
xk_lqr = xerr(1,:)'*0;
for kk = 1:N*2
%     yvec(kk) = sys_sim.c*xk;
    xvec_est(:,kk) = xk;
    uk = K_est*(xf - xk);
    uk = max(min(uk, du_max), -du_max);
    uvec_est(kk) = uk;
    xk = sys_sim.a*xk + sys_sim.b*uk;

    % -------------- Now LQR ---------------%
    xvec_lqr(:,kk) = xk_lqr;
    uk_lqr = K_lqr*(xf - xk_lqr);
    % (norm(K_est)/norm(K_lqr))*
    uk_lqr = max(min(uk_lqr, du_max), -du_max);
    uvec_lqr(kk) = uk_lqr;
    xk_lqr = sys_sim.a*xk_lqr + sys_sim.b*uk_lqr;

    
    
end


figure(1);
subplot(2,1,1)
hold on
plot(sys_sim.c*xvec_est, '--')
plot(sys_sim.c*xvec_lqr, ':')


subplot(2,1,2)
hold on
uvec_est = max(min(uvec_est, du_max),-du_max);
plot(uvec_est, '--')

plot(uvec_lqr, ':')



% plot(k_usat, u_sat, 'xk')
%%
kkdot = -K_lqr*K_est'
costheta = norm(K_lqr, 2)*norm(K_est, 2)/kkdot
acosd(costheta)



