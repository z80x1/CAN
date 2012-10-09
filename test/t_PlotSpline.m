function PlotSpline % ������ ������ �������
global mip scw sch x y xx hFig hAxes;
yy = spline(x, y, xx);
% �������� ������� �������
axes(hAxes); % ���������� ������� ���
cla % �������� ���
hold on % �������� ��������
hPoint = line(x, y); % ������� �����
hLine = line(xx, yy); % ������
set(hPoint, 'LineStyle', 'none', 'Marker',...
'square', 'MarkerSize', 7, 'MarkerEdgeColor',...
[0 0 1], 'MarkerFaceColor', [0 0 1]);
set(hLine, 'Color', [0 0 0]);
hold off % ��������� ��������
return % ������� � ���������� ���������