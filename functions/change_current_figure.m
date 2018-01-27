% change_current_figure(h)
% 
% Change current figure without changing focuse to it;
%
% Inputs
% ------
%   h :  Handle to the figure you want to plot to.

function change_current_figure(h)
    set(0,'CurrentFigure',h)
end
