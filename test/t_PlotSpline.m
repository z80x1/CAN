function PlotSpline % рисует график сплайна
global mip scw sch x y xx hFig hAxes;
yy = spline(x, y, xx);
% ординаты графика сплайна
axes(hAxes); % установили текущие оси
cla % очистили оси
hold on % включили задержку
hPoint = line(x, y); % узловые точки
hLine = line(xx, yy); % сплайн
set(hPoint, 'LineStyle', 'none', 'Marker',...
'square', 'MarkerSize', 7, 'MarkerEdgeColor',...
[0 0 1], 'MarkerFaceColor', [0 0 1]);
set(hLine, 'Color', [0 0 0]);
hold off % выключили задержку
return % возврат в вызывающую программу