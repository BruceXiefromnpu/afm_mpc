
% Ts = 0.01;
% sys = c2d(ss(tf([144], [1, 2*12*0.01, 144])), Ts);
% sys.InputDelay =0;
% sys = absorbDelay(sys);
% % figure(5)
% % bode(sys)
% %
% wnyq = (pi/Ts);
% omegas1 = logspace(log10(0.005), log10(wnyq), 200);
% [~,~, omegas] = bode(sys);
% omegas = unique([omegas1(:); omegas(:)] );
% 
% 
% 
% sys_recyc = SSTools.deltaUkSys(sys);
% pint = 0.8;
% Px = getCharDes(sys_recyc, [1, 1], pint, [0.85, 0.85], [2], .25);
% [chat, dhat] = place_zeros(sys_recyc, Px);
% Q = chat'*chat;
% S = chat'*dhat;
% R0 = dhat^2;
% 
% 
% 
% 
% a = sys_recyc.a;
% b = sys_recyc.b;
% [Ns, Nu] = size(b);


function K = place_varga(a, b Px)
clear classes
mod = py.importlib.import_module('call_pc');
py.reload(mod);

if count(py.sys.path, '') == 0
  insert(py.sys.path, int32(0), '')
end

[Ns, Nu] = size(b);

K = py.call_pc.place_varg(a(:)', b(:)', int32(Ns), int32(Nu), real(Px), imag(Px));
K = double(K);


end