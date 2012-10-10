function my_axis(limits)
%scale axis with 1% borders and allows NaN values
min_x = limits(1);
max_x = limits(2);
min_y = limits(3);
max_y = limits(4);

        ysize=max_y-min_y;
        if isnan(ysize)
            min_y = 0;
            max_y = 1;
        else
            min_y = min_y-0.01*ysize;
            max_y = max_y+0.01*ysize;
        end
        xsize=max_x-min_x;
            min_x = min_x-0.01*xsize;
            max_x = max_x+0.01*xsize;

        axis([min_x max_x min_y max_y]);
                                           