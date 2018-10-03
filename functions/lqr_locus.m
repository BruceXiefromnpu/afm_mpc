

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
  
  hands = gobjects(length(gam_s), 1);
  for k=1:length(gam_s)
    R = Ro + gam_s(k);
    
    K_lqr = dlqr(sys.a, sys.b, Q, R, S);
    sys_cl = SSTools.close_loop(sys, K_lqr);
    
    pl = pole(sys_cl);
    clr = mp_fine(k,:);
    hands(k) = plot(real(pl), imag(pl), '.', 'color', clr);
    hands(k).UserData = gam_s(k);
    hold on
    drawnow();
  end
  dcm = datacursormode(gcf);
  set(dcm, 'UpdateFcn', @local_callback_function);
    
%   tcks = logspace(log10(gam_s(1)), log10(gam_s(end)), 10);
%   tcks =  round(tcks, 2, 'significant');
% 
%   caxis([gam_s(1), gam_s(end)])
%   C_hand = colorbar('TickLabels', tcks, 'Location', 'northoutside'); 
  
  xlabel('Re');
  ylabel('Im');
  
end



function datatip_txt = make_tip_txt(sz_loc, freq, Wn, zet, mag)
if abs(Wn - freq) < 1e-13
    unit = '[rad/s]';
else
    unit = '[Hz]';
end
    datatip_txt = {['Freq: ',num2str(freq,4), ' ',unit],...
        ['Mag: ',num2str(mag,4), ' [dB]'],...
        ['Wn: ',num2str(Wn,4), ' rad/s'],...
        ['damping: ',num2str(zet,4)],...
        ['Re: ', num2str(real(sz_loc))],...
        ['Im: +-', num2str(abs(imag(sz_loc)))]};
end

function output_txt = local_callback_function(~,event_obj)
% Display the position of the data cursor
% obj          Currently not used (empty)
% event_obj    Handle to event object

    if isa(event_obj.Target.UserData, 'double')
      % check to see if the target line has the user data set. If so,
      % assume it will return the text to update the data-tip.
      gam = event_obj.Target.UserData();
    else
      gam = NaN;
    end
    
    pos = get(event_obj, 'Position');
    
    output_txt = {['X: ',num2str(pos(1))],...
      ['Y: ',num2str(pos(2))],...
      ['gam: ', num2str(gam)]};


end




