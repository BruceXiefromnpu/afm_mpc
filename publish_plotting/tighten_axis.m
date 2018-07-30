function tighten_axis(fig, ax)
%   tighten_axis(fig, ax)
% 
% This function only works for a single axis. It will try to pull
% the axis limits as tight as possible to the edges of the figure.
% 
% 
% The legend size and position is a bit weird because matlab
% doesn't always honor what you set the Position property to. So
% for now, this function does not touch the legend. Thus, in a
% typical workflow, you would build the figure without messing with
% the axis sizes etc, call
% tighten_axis(fig, ax)
% then insert your legend and set its position exactly where you
% want it.
% 
% 
% Notes:
% There are several coordiniate systems. 
% Figure:
%   (PaperPosition). Since we are doing this to make publication
%   figures, ensure that Position and PaperPosition are the same,
%   at least in width and height and units.
%   Position
% Axis:
%   Position. This is relative to the Figure.Position. 
% Text Positions, like title, xlabel, ylabel. These are relative to
% the axis.This means that the vertical position (y-coord) of the
% xlabel is typically negative and the horizontal position of the
% y-label is typically negative. 
% 
% This function assumes that the xticks and yticks extend to the
% edges of the axis. Generally, the horizontal middle of the xtick
% is aligned with the left axis edge, and the vertical middle of
% the top ytick is aligned with top axis edge. This function tries
% to account for that measuring the size of the last xtick and
% ytick and spacing the axis left (down) by half the width (height)
% of that text. 
  set(fig, 'Units', 'inches')
  set(ax, 'Units', 'inches', 'FontUnits', 'points')
  
  pos = fig.Position;
  paper_pos = fig.PaperPosition;
  if pos(3) ~= paper_pos(3) || pos(4) ~= paper_pos(4)
    error(['Please ensure paperPosition and Position are equal. ' ...
           'Otherwise incostient, confusing results will ensure. ' ...
           'Consider using mkfig(num, width, height) for this'])
  else
  
    fig_height = paper_pos(4);
    fig_width  = paper_pos(3);
    
  end
  
  set(findobj(fig, 'Tag', 'legend'), 'Units', 'inches', 'FontUnits', 'points');
                    
  yt_pad = 1.3; % how much we imagine the yticks are spaced from
                % the y-axis.
  
  ylab = ax.YLabel;
  xlab = ax.XLabel;
  tit_lab = ax.Title;
  
  set(xlab, 'Units', 'inches', 'FontUnits', 'points');
  set(ylab, 'Units', 'inches', 'FontUnits', 'points');
  set(tit_lab, 'Units', 'inches', 'FontUnits', 'points');
  
  
  [ylab_wd, ylab_ht] = get_font_wd_ht(ylab);
  [xlab_wd, xlab_ht,] = get_font_wd_ht(xlab);
  [~, tit_lab_ht] = get_font_wd_ht(tit_lab);
  
  [xt_wd, xt_ht] = get_xticks_wd_ht(ax);
  [yt_wd, yt_ht] = get_yticks_wd_ht(ax);
  % compute the bottom of the axes. Account for height of xlabel,
  % height of xticks, plus 30% for space between xticks and axis
  
  bt = xlab.FontSize/72 + xt_ht;
  lft = ylab_wd + yt_pad*yt_wd;
  
  % compute axis height. This is the height of the figure, minus
  % room on the bottome (bt) for xlabel and xticks, minus room on
  % the top for title and/or half of the ytick height (since the
  % top ytick is centered around the top of the axis box).
  
  ht = fig_height - bt - max(tit_lab_ht, 0.5*yt_ht);
  
  wd = fig_width - lft - xt_wd/2;
  
  set(ax, 'Position', [lft, bt, wd, ht])
  set(xlab, 'Position', [ax.Position(3)/2, -xlab.FontSize/72, 0])

  set(ylab, 'HorizontalAlignment', 'center', 'Position', [-yt_wd*yt_pad, ax.Position(4)/2, 0])
  

end


function [width, height] = get_xticks_wd_ht(ax)
  xtick_max_str_s = get(ax, 'XTickLabels');
  
  FS = get(ax, 'FontSize');
  
  % There doesn't seem to be any access to the yticks and xticks.
  % Thus, create a text object with the largest string, measure it,
  % then delete it. 
  
  h_temp = text(1,1, xtick_max_str_s{end}, 'FontSize', FS, 'Units', 'inches');
  
  width = h_temp.Extent(3);
  height = h_temp.Extent(4);
  
  delete(h_temp)
  
end


function [width, height] = get_yticks_wd_ht(ax)

  % Take max by string length, not number, because '10' takes up less space
  % that '0.001'
  ytick_str_s = get(ax, 'YTickLabels');
  nticks = length(ytick_str_s);
  lens = zeros(1, nticks);
  for i=1:nticks
    lens(i) = length(ytick_str_s{i});
  end
  [~, ytick_max_idx] = max(lens);
  FS = get(ax, 'FontSize');
  
  % There doesn't seem to be any access to the yticks and xticks.
  % Thus, create a text object with the largest string, measure it,
  % then delete it. 
  h_temp = text(1,1, ytick_str_s{ytick_max_idx}, 'FontSize', FS, 'Units', 'inches');
    
  width = h_temp.Extent(3);
  height = h_temp.Extent(4);
  
  delete(h_temp)
  
end


function [width, height] = get_font_wd_ht(fig_str)
% Returns the height and width taken by  figure-string object in inches.
  units = get(fig_str, 'FontUnits');
  set(fig_str, 'FontUnits', 'points', 'Units', 'inches');
  
  ext = fig_str.Extent();
  height = ext(4);
  width = ext(3);

  set(fig_str, 'FontUnits', units)
  
end
