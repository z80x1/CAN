function my_axis(params, plottype, cutoffs, inlineindex, what_rescale, fl_plot_numbers)
%scale axis with 1% borders and allows NaN values

    global flags

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
        inlineindex=1;
    end
    if nargin<5
        what_rescale=[1 1 1 1];
    end
    if nargin<6
        fl_plot_numbers=0;
    end

    limits = params.limits;
    min_x = limits(1);
    max_x = limits(2);
    min_y = limits(3);
    max_y = limits(4);
    borders=[0.05 0.05 0.02 0.10];

    fontsize = 11;
    axis_fontsize = 15;
    axislabel_fontsize = 17;

%limits
    ysize=max_y-min_y;
    if what_rescale(3)
        if isnan(ysize)
            min_y = 0;
        else
            min_y = min_y-borders(3)*ysize;
        end
    end
    if what_rescale(4)
        if isnan(ysize)
            max_y = 1;
        else
            max_y = max_y+borders(4)*ysize;
        end
    end
    
    xsize=max_x-min_x;
    if what_rescale(1)
        min_x = min_x-borders(1)*xsize;
    end
    if what_rescale(2)
        max_x = max_x+borders(2)*xsize;
    end
    axis([min_x max_x min_y max_y]);

%cutouffs    
    cutoffs=cutoffs(~isnan(cutoffs));
    cutoffs=sort(cutoffs);
    if numel(cutoffs)
        plot(repmat(cutoffs,2,1), repmat([min_y; max_y],1,numel(cutoffs)), 'k', 'LineWidth', 0.25); %plot cutoff lines
        for i=1:numel(cutoffs)
            ypos = max_y+(-0.07+0.045*mod(i,2))*ysize;
            ht = text(cutoffs(i), ypos, int2str(i), 'FontSize',fontsize, 'Rotation', 0, 'BackgroundColor','w', 'Margin',eps);
%                text(cutoffs(i), ypos, int2str(i), 'FontSize',6, 'Rotation', 90, 'BackgroundColor','w' )
            if fl_plot_numbers && flags.develmode
%                ypos = min_y+(0.18+0.2*mod(i,2))*ysize;
                ypos = max_y+(-0.20-0.15*mod(i,2))*ysize;
                ht = text(cutoffs(i), ypos, num2str(cutoffs(i),'%0.2f'), 'FontSize',fontsize, 'Rotation', 90, 'BackgroundColor','w', 'Margin',eps);
            end
        end
    end


    %ticks
    if ~params.flags.develmode && ~isnan(params.flags.irc_tick_interval)
        xticks = ceil(min_x):params.flags.irc_tick_interval:floor(max_x);
        set(gca,'XTick', xticks);
    else
        set(gca,'XTickMode','auto');
    end
    
    if strncmp(plottype,'rho',3) || strncmp(plottype,'deltarho',8);
        ticks_format('','%0.3f');
    elseif strncmp(plottype,'epsilon',7);
        ticks_format('','%0.2f');
    elseif strcmp(plottype,'dHH');
        ticks_format('','%0.1f');
    else
        set(gca,'YTickLabelMode','auto');
    end
    
    set(gca,'Box','on','XMinorTick','on','YMinorTick','on','FontSize',axis_fontsize);

    %labels
    xlabel_str = 'IRC, Bohr';
    hx = xlabel(xlabel_str,'interpreter','latex','FontSize',axislabel_fontsize); %#ok<NASGU>

    interpreter='latex';
    switch plottype 
        case 'E'
            ylabel_str='$\Delta$ E, kcal/mol'; 
        case 'diffE_irc'
            ylabel_str='{dE}/{dIRC}';
        case 'mu'
            ylabel_str='$\mu$, D';
        case {'dAB','dAB_CHB'}
            ylabel_str='$d_{AB}, \AA$';
        case {'dHB'}
            ylabel_str='$d_{AH/HB}, \AA$';
        case {'dHB_CHB'}
            ylabel_str='$d_{HB}, \AA$';
        case {'aAHB','aAHB_CHB'}
            ylabel_str='{\angle}AH...B, degree';
%             labelStr = '<html>&#x2220;</html>'; % angle symbol
%             jLabel = javaObjectEDT('javax.swing.JLabel',labelStr);
%             [hcomponent,hcontainer] = javacomponent(jLabel,[100,100,20,20],gcf);
        case 'Milliken'
            ylabel_str='q, e';
        case 'NBO_q'
            ylabel_str='NBOq_{H}, e';
        case 'dHH'
            ylabel_str='R(H-H), \AA';
        case 'glycang'
            ylabel_str='$\alpha$, degree';
        case {'rho','rhoCHB'}
            ylabel_str='$\rho$, a.u.';
        case {'deltarho','deltarhoCHB'}
            ylabel_str='$\Delta \rho$, a.u.';
        case {'epsilon','epsilonCHB'}
            ylabel_str='$\epsilon$';
        case {'EHB','EHBCHB'}
            ylabel_str='E_{HB}, kcal/mol';interpreter='tex';
        case 'EHB_timn'
            ylabel_str='E_{HB} Nikolayenko, kcal/mol';
        otherwise
            ylabel_str='undefined';
    end
    hy = ylabel(ylabel_str,'interpreter',interpreter,'FontSize',axislabel_fontsize);
    if inlineindex~=1
        set(hy,'Color',[1 1 1]);
    end
    
    ht = get(gca,'title');
    delete(ht);
    grid off

                                                                