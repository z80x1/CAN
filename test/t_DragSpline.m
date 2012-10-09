clear all % очистили рабочую область
global mip scw sch x y xx hFig hAxes;

scrsz = get(0, 'ScreenSize'); % размеры экрана монитора
wscr = 640; % ширина окна
hscr = 440; % высота окна
waxe = wscr-40; % ширина осей
haxe = hscr-40; % высота осей

wgr = 6; % диапазон по x
hgr = 4; % диапазон по y

scw = waxe/wgr; % масштаб по x
sch = haxe/hgr; % масштаб по y
pl = (scrsz(3)-wscr)/2; % левый край фигуры
pu = (scrsz(4)-hscr)/2; % верхний край фигуры
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

x = 0:wgr; % абсциссы узловых точек
y = rand(1, wgr+1)*hgr; % ординаты узловых точек
xx = 0:0.01:wgr; % абсциссы точек графика сплайна

t_PlotSpline; % рисуем график сплайна