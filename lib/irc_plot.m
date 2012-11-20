function h = irc_plot( AX, x, y, initial_index )
%AX - axis to plot to
%x, y - data to plot
%initial_index - plot number on the axis

persistent currentCount;
global flags

if nargin == 4
    currentCount = initial_index;
else
    currentCount = currentCount + 1;
end


pcolor=[{'r'} {'g'} {'b'} {'k'} {'m'} {'c'} ...
        {[.5412 .1686 .8863]}  ...
        {[0 .5 0]} {[.5 0 0]} {[0 0 .5]} {[.5 .5 0]} {[.5 0 .5]} {[0 .5 .5]} {[.5 .5 .5]} ...
        {[.5 .1 .9]} {[.8 .2 .2]} {[.8 .8 .2]}...
        {[.9 .4 .9]} {[.2 .4 .6]} {[.6 .4 .6]} {[.6 .2 .2]} {[.8 .2 .8]} ...
        {[.2 .8 .8]} {[.2 .8 .2]} {[.2 .2 .8]} {[.4 .9 .1]} {[.1 .3 .6]} {'y'} {[.75 .75 .75]} {[.2745 .5098 .7059]}];
psign='x.ov^x<>d*hp+';
%psign='ov^x<>d*hp+.';
pstyle=[{'-'} {'--'} ];
lw = 0.5; %line width
ms = 3; %marker size (points)
mscrosses = 5; %marker size for crosses (points)
ms_dots = 5; %marker size for dots


pointsign = psign(mod(currentCount-1,numel(psign))+1);
if pointsign=='x' %crosses are too small - lets increase their size 
    markersize = mscrosses;
elseif pointsign=='.' %dots are too small - lets increase their size 
    markersize = ms_dots;
else
    markersize = ms;
end
pointsign = [pointsign '-'];

if flags.develmode
    color = pcolor{mod(currentCount-1,numel(pcolor))+1};
else    
    color = 'k';
end

if ~flags.develmode
    %find gap in plot and filling it
    avg_dx = (max(x)-min(x))/numel(x);
    min_gap_width = 5*avg_dx;

    dx = x(2:end)-x(1:end-1);
    igap = find(dx>min_gap_width);
    if ~isempty(igap)
        disp(['gap in plot found. filling it']);
        xgap = [x(igap) x(igap+1)];
        ygap = [y(igap) y(igap+1)];
        xg = xgap(1):avg_dx:xgap(2);
        yg = (xg-xgap(1)) .* (ygap(2)-ygap(1))/(xgap(2)-xgap(1)) + ygap(1);

        [x,isort] = sort([x xg]);
        y = [y yg];
        y = y(isort);
    end    
end

h = plot( AX, x,y,pointsign,'Color',color,'MarkerFaceColor',color,'LineWidth',lw,'MarkerSize', markersize);

