Vdiv = 2.56/(19.61+2.56);
Vdiv_gain = 1/Vdiv;

load( fullfile(PATHS.exp, 'x-axis_sines_info_out_2-8-2018-01.mat'))
Gpow = ss(modelFit.models.G_uz2pow*Vdiv_gain);
Gpow.InputDelay = 0


%%

clc
dumax = 1;
k0 = 78;

PHI = sys_sim.a;
Gam = sys_sim.b;
PHI_pow = Gpow.a;
Gam_pow = Gpow.b;
C_pow = Gpow.c;

H=[];
H_pow = []

F = [];
phi = [];
Ns = 2; nu = 1;
no=1;

for i = 1:k0
    phi = [phi; C_pow*PHI_pow^i];
end

F_1 = [C_pow*eye(Ns); phi(1:end-no, :)]*Gam_pow;

for i=0:k0-1
   F = [F, [zeros(no*i,nu); F_1(1:no*k0-(i)*no, :)] ];
end
%
for k=0:k0-1
    H = [PHI^k*Gam H];
%     H_pow = [C_pow*PHI_pow^k*Gam_pow H_pow];
end
e_k0 = zeros(1,k0);
e_k0(end) = 1; % Kth unit vector. [0 0 .... 0 0 1]. To pick out last element of u.
I = eye(length(Gam));

Fcumsum = zeros(k0,k0);
for k=1:k0
   Fcumsum(k, 1:k) = ones(1,k); 
end

% variable u(k0);

% L = self.PHI^k0*x0 + C*u;
opts = optimset(@lsqlin);
opts.Algorithm = 'interior-point';
% opts.Algorithm = 'active-set';
opts.MaxIter = 20000;
d = -PHI^k0*x0 + xss_sim;
ub = ones(k0,1)*dumax;
[u2,~,res] = lsqlin(H, d, [F;-F; Fcumsum;-Fcumsum],[ub;ub; 0*ub+7; 0*ub+7], [],[], [],[] ,[], opts);
res
if isvalid(h2)
    delete(h2)
end

figure(2),clf, hold on, grid on;
h2 = plot(u2, 'r');

uu2 = [u2; zeros(400,1)];
t2 = [0:1:length(uu2)-1]'*Ts;
uu = lsim(g_eject, uu2, t2);
%
y = lsim(SYS, cumsum(uu), t2);
% y = lsim(sys_sim, uuu, [0:1:k0-1]'*Ts);

figure(1);clf, hold on
plot(t2, cumsum(uu), '--g')
plot(t2, y, '-r')

ypow = lsim(Gpow, uu2, t2);
figure(10);
plot(t2, ypow)

