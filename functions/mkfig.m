function F = mkfig(num, width, height)

F = figure(num); clf
axes('FontName','Times New Roman') % Set axis font style
box('on'); % Define box around whole figure
hold on;

set(F, 'Units', 'Inches', 'Position', [0, 0, width, height],...
    'PaperUnits', 'Inches', 'PaperSize', [width, height])
set(F, 'Color', 'w');



end