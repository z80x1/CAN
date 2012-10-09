clear all % �������� ������� �������
global mip scw sch x y xx hFig hAxes;

scrsz = get(0, 'ScreenSize'); % ������� ������ ��������
wscr = 640; % ������ ����
hscr = 440; % ������ ����
waxe = wscr-40; % ������ ����
haxe = hscr-40; % ������ ����

wgr = 6; % �������� �� x
hgr = 4; % �������� �� y

scw = waxe/wgr; % ������� �� x
sch = haxe/hgr; % ������� �� y
pl = (scrsz(3)-wscr)/2; % ����� ���� ������
pu = (scrsz(4)-hscr)/2; % ������� ���� ������
hFig = figure('Position',[pl pu wscr hscr],...
'Resize', 'off', 'MenuBar', 'none',...
'NumberTitle', 'off',...
'WindowButtonMotionFcn', 't_ChangePointer;',...
'WindowButtonDownFcn',  'mip=1;',...
'WindowButtonUpFcn',  'mip=0;',...
'Name', 'Spline Interpolation');

hAxes = axes('Units', 'pixels', 'Position',...
[20 20 waxe haxe], 'XLim', [0 wgr], 'YLim',...
[0 hgr]);

x = 0:wgr; % �������� ������� �����
y = rand(1, wgr+1)*hgr; % �������� ������� �����
xx = 0:0.01:wgr; % �������� ����� ������� �������

t_PlotSpline; % ������ ������ �������