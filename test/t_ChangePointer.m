function t_ChangePointer % ������ � ������
global mip scw sch x y hFig ;

p = get(hFig, 'CurrentPoint'); % ���������� �����
xscr = x*scw+20; % �������� �������� ������� �����
yscr = y*sch+20; % �������� �������� ������� �����
pinp = find((p(1)>xscr-3)&(p(1)<xscr+3)&...
    (p(2)>yscr-3)&(p(2)<yscr+3)); % ����� ���� ���������?

if isempty(pinp), % ��� ����������� �����
 set(hFig, 'Pointer', 'arrow'); % ���������  �������
else % ���� ����������� ����
 set(hFig, 'Pointer', 'crosshair'); % ���������  �������
end

if (~isempty(pinp))&mip, % ������ ������ � ����� ���������
  x(pinp) = (p(1)-20)/scw; % ����� �������� ����
  y(pinp) = (p(2)-20)/sch; % ����� �������� �������������� � �����������
  t_PlotSpline;   % ������������ ������ � ������� �����
end
return % ������� � ���������� ���������