function h = irc_plot( AX, x, y, initial_index )

persistent currentCount;

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
    psign='ov^x<>d*hp+.';
    pstyle=[{'-'} {'--'} ];
    lw = 0.5; %line width
    ms = 2; %marker size (points)
    mscrosses = 3; %marker size for crosses (points)


                 pointsign = psign(mod(currentCount-1,numel(psign))+1);
                 if pointsign=='x' %crosses are too small - lets increase their size 
                     markersize = mscrosses;
                 else
                     markersize = ms;
                 end
                 pointsign = [pointsign '-'];

color = pcolor{mod(currentCount-1,numel(pcolor))+1};


h = plot( AX, x,y,pointsign,'Color',color,'MarkerFaceColor',color,'LineWidth',lw,'MarkerSize', markersize);    

