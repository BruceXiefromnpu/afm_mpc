

function [ax, C_hand] = lqr_locus(sys, Q, Ro, S, lb, ub, ax, varargin)
  
  
  if ~exist('ax', 'var') || isempty(ax)
    ax = gca;
  elseif ~isa(ax, 'matlab.graphics.axis.Axes')
    error('ax must be a graphics axis\n')
  end
  
  p = inputParser;
%   p.addParameter('ax', gca);
  p.addParameter('color_map', []);
  p.addParameter('N', 100);
  p.parse(varargin{:});
  N = p.Results.N;
  
  mp = p.Results.color_map;
  if isempty(mp)
    mp = colormap(ax, 'jet');
  else
    colormap(ax, mp)
  end
%   ax = p.Results.ax;
  
  
  if isempty(S)
    S = sys.b*0;
  end
  
  gam_s = logspace(log10(lb), log10(ub), N);

  x = 1:size(mp,1);

  xq = linspace(1, (size(mp,1)), length(gam_s));

  mp_fine = zeros(length(xq), 3);

  mp_fine(:,1) = interp1(x, mp(:,1), xq);
  mp_fine(:,2) = interp1(x, mp(:,2), xq);
  mp_fine(:,3) = interp1(x, mp(:,3), xq);
  
  pl = pole(sys);
  plot(real(pl), imag(pl), 'xk', 'MarkerSize', 10)
  hold on;
  if sys.Ts > 0
    zgrid();
  else
    sgrid();
  end
  
  for k=1:length(gam_s)
    R = Ro + gam_s(k);
    
    K_lqr = dlqr(sys.a, sys.b, Q, R, S);
    sys_cl = SSTools.close_loop(sys, K_lqr);
    
    pl = pole(sys_cl);
    clr = mp_fine(k,:);
    plot(real(pl), imag(pl), '.', 'color', clr)
    
    hold on
    drawnow();
  end
  
  tcks = logspace(log10(gam_s(1)), log10(gam_s(end)), 10);
  tcks =  round(tcks, 2, 'significant');

  caxis([gam_s(1), gam_s(end)])
  C_hand = colorbar('TickLabels', tcks, 'Location', 'northoutside'); 
  
  xlabel('Re');
  ylabel('Im');
  
end