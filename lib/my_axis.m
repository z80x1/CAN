function my_axis(limits, plottype, cutoffs, what_rescale, fl_plot_numbers)
%scale axis with 1% borders and allows NaN values

if nargin<1
    error('no arguments');
end
if nargin<2
    plottype='auto';
end
if nargin<3
    cutoffs=[];
end
if nargin<4
    what_rescale=[1 1 1 1];
end
if nargin<5
    fl_plot_numbers=0;
end
cutoffs=cutoffs(find(~isnan(cutoffs)));
cutoffs=sort(cutoffs);

    min_x = limits(1);
    max_x = limits(2);
    min_y = limits(3);
    max_y = limits(4);
    border=0.05;

    ysize=max_y-min_y;
    if what_rescale(3)
        if isnan(ysize)
            min_y = 0;
        else
            min_y = min_y-border*ysize;
        end
    end
    if what_rescale(4)
        if isnan(ysize)
            max_y = 1;
        else
            max_y = max_y+border*ysize;
        end
    end
    xsize=max_x-min_x;
    if what_rescale(1)
        min_x = min_x-border*xsize;
    end
    if what_rescale(2)
        max_x = max_x+border*xsize;
    end

    if numel(cutoffs)
        plot(repmat(cutoffs,2,1), repmat([min_y; max_y],1,numel(cutoffs)), 'k')
            for i=1:numel(cutoffs)
                ypos = min_y+(0.1+0.2*mod(i,2))*ysize;
                text(cutoffs(i), ypos, int2str(i), 'FontSize',6, 'Rotation', 90, 'BackgroundColor','w' )
                if fl_plot_numbers
                    ypos = min_y+(0.18+0.2*mod(i,2))*ysize;
                    text(cutoffs(i), ypos, num2str(cutoffs(i),'%0.2f'), 'FontSize',6, 'Rotation', 90, 'BackgroundColor','w' )
                end
            end
    end

            
    axis([min_x max_x min_y max_y]);
    if strncmp(plottype,'rho',3);
        ticks_format('','%0.3f');
    elseif strncmp(plottype,'deltarho',8);
        ticks_format('','%0.3f');
    elseif strncmp(plottype,'epsilon',7);
        ticks_format('','%0.2f');
    else
        set(gca,'YTickLabelMode','auto');
    end
    
                                                              