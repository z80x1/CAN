function ChangePointer(hFig,mol) % работа с мышкой
global mip hFig;

disp('Hi!');
p = get(hFig, 'CurrentPoint'); % координаты мышки

xscr = mol.x; % экранные абсциссы узловых точек
yscr = mol.y; % экранные ординаты узловых точек
pinp = find((p(1)>xscr-3)&(p(1)<xscr+3)&...
    (p(2)>yscr-3)&(p(2)<yscr+3)); % какие узлы захвачены?

if isempty(pinp), % нет захваченных узлов
 set(hFig, 'Pointer', 'arrow'); % указатель  стрелка
else % есть захваченный узел
 set(hFig, 'Pointer', 'crosshair'); % указатель  крестик
end

%if (~isempty(pinp))&mip, % нажата кнопка и точка захвачена
%  x(pinp) = (p(1)-20)/scw; % нова€ абсцисса узла
% y(pinp) = (p(2)-20)/sch; % нова€ ордината узлаћатематика в приложени€х
% t_PlotSpline;   % перерисовали сплайн и узловые точки
%nd
return % возврат в вызывающую программу