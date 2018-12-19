clear



g = zpk([0.25], [.5+j, .5-j, -20], 22);
[gm0, pm0, wcg0, wcp0] = margin(g);
[re, im] = nyquist(g);
figure(4); clf
ax1 = subplot(2,2,1)
plot(re(:), im(:))
hold on, grid on
plot(re(:), -im(:), '--')
plot_nyq_pm(ax1, g)
legend('$\omega>0$', '$\omega <0$')
title(sprintf('G(jw): GM = %.2f, PM $= %.1f^{\\circ}$ (CL stable)', gm0, pm0))

subplot(2,2,[3,4])
bode(g)
grid on, hold on

% title(sprintf('Is stable?=%d', isstable(feedback(g, 1))))
%
% Design a derivative controller which has a gain of 1 at wcp, and provides
% more than 15.5 degrees of phase.
phi = 17 *pi/180;
alp = wcp0 / tan(phi);

% compute gain beta such that beta *(jwo + alp) = 1
beta = 1/sqrt(wcp0^2 + alp^2);

K = zpk(-alp, [], beta);
abs(freqresp(K, wcp0))

L = g*K
[gm1, pm1, wcg1, wcp1] = margin(L);

[re, im] = nyquist(L);
re = re(:);
im = im(:);

ax2 = subplot(2,2,2);
plot(re, im)
hold on
plot(re, -im, '--')
plot_nyq_pm(ax2, L)
title(sprintf('D(jw)G(jw): GM = %.2f, PM $= %.1f^{\\circ}$ (CL unstable)', gm1, pm1))
grid on

subplot(2,2,[3,4])
hold on
bode(g, L)
legend('$G(s)$', '$D(s)G(s)$')
grid  on
axs = get(gcf, 'Children');
%%
ylim(axs(4), [-80, 5])


function plot_nyq_pm(ax, L)
[gm, pm, wcg, wcp] = margin(L);
gm_abs = 10^(gm/20);

% idx_wcg = w(w==wcg);
% idx_wcp = w(w==wcp);
L_wcp = squeeze(freqresp(L, wcp));
L_wcg = squeeze(freqresp(L, wcg));

t = [0:0.01:2*pi]';
x = cos(t);
y = sin(t);
plot(x, y, ':')
plot(ax, [-1], 0, 'rx', 'MarkerSize', 5);
plot(ax, [0 real(L_wcp)], [0, imag(L_wcp)], '-k')
plot(ax, [real(L_wcp)], [imag(L_wcp)], '.', 'MarkerSize', 15)

try
  H = feedback(L, 1);
  flg = isstable(H);
catch
  figure(100)
  [y, t] = step(H);
  plot(t, y);
  figure(ax.Parent);
  if max(abs(y)) < 1000
    flg = true;
  else
    flg = false;
  end
end

if flg
  stab_str = 'yes';
else
  stab_str = 'no';
end

fprintf('GM = %f.   PM = %f\n Closed-Loop stable = %s\n', gm_abs, pm, stab_str)

end