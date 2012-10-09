function ChangePointer(hFig,mol) % ������ � ������
global mip hFig;

disp('Hi!');
p = get(hFig, 'CurrentPoint'); % ���������� �����

xscr = mol.x; % �������� �������� ������� �����
yscr = mol.y; % �������� �������� ������� �����
pinp = find((p(1)>xscr-3)&(p(1)<xscr+3)&...
    (p(2)>yscr-3)&(p(2)<yscr+3)); % ����� ���� ���������?

if isempty(pinp), % ��� ����������� �����
 set(hFig, 'Pointer', 'arrow'); % ���������  �������
else % ���� ����������� ����
 set(hFig, 'Pointer', 'crosshair'); % ���������  �������
end

%if (~isempty(pinp))&mip, % ������ ������ � ����� ���������
%  x(pinp) = (p(1)-20)/scw; % ����� �������� ����
% y(pinp) = (p(2)-20)/sch; % ����� �������� �������������� � �����������
% t_PlotSpline;   % ������������ ������ � ������� �����
%nd
return % ������� � ���������� ���������