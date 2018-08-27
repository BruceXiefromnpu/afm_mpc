
clear
clc

g = load('gdrift_inv.mat');
g = canon(g.gdrift_inv);

% Normalize g.b to one.
nx = SSTools.getNxNu(g);

T = diag(1./nx);

g2 = ss2ss(g, T);


t = (0:floor(0.04/g2.Ts))'*g2.Ts;
u = t*0 + 1;
[y, t, x] = lsim(g2, u, t);
nw = 18;
a = sfi(g2.a, nw, nw-1);
b = sfi(g2.b, nw, nw-1);
c = sfi(g2.c, nw, nw-2);
d = sfi(g2.d, nw, nw-2);


xk = sfi([0;0], nw, nw-5);

uu = sfi(u, nw, nw-5);
yy = uu*0;
xfxp = [];
for k=1:length(uu)

  yy(k) = c*xk + d*uu(k);
  
  xk = sfi(a*xk + b*uu(k), nw, nw-5);
  
  xfxp = [xfxp; xk'];
    
end



figure(1); clf
plot(t, x)
hold on
plot(t, xfxp, '--')

figure(2); clf;
plot(t, y, '-b', t, yy, '--r')