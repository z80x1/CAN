function res = movetexts(buf,flags)
%try to organize nonoverlapping plot labels 
res=0;

if ~numel(buf.ht)
    return % nothing to do
end
fl_debug = 0;

exts={}; y=[];
for it=1:numel(buf.ht)
    exts(it).ht = buf.ht(it);
    exts(it).h_plots = buf.h_plots(it);
    exts(it).extent = get(buf.ht(it),'extent');
    exts(it).success = 0;
    exts(it).str = get(buf.ht(it),'string');
    exts(it).fontsz = get(buf.ht(it),'FontSize');
    exts(it).disabled = 0;
end
%childs=get(gca,'children')
x = get(buf.h_plots(1),'XData');
for ip=1:numel(buf.h_plots)
    y(ip,:) = get(buf.h_plots(ip),'YData');
    
    maxy = max(abs(y(ip,:)));
    miny = min(abs(y(ip,:)));
    for ip2=1:ip-1
        maxdelta = max(abs(y(ip2,:)-y(ip,:)));
        
        if (maxy/maxdelta)>39 %similar plots %39 is critical value for Hyp-Hyp E_HB O6H
            exts(ip).disabled = 1; %ip is same as it
        end
    end
end

ixstep = floor(size(x,2)/20); %analyze every 20th point
ipoints = 2:ixstep:numel(x); %start to deal from second point 
[dummy,i_from_ends2center]=sort(ipoints.*(max(ipoints)-ipoints)); %indexes for moving from ends to center

xlimits=xlim;
ylimits=ylim;
ymgn_ownline = 0.01*diff(ylimits);
xmgn_alien = 0.01*diff(xlimits);
ymgn_alien = 0.02*diff(ylimits);
xmgn_box = 0.05*diff(xlimits);
ymgn_box = 0.02*diff(ylimits);

garbobjs = [];

for it=1:numel(exts)
    if exts(it).disabled == 1, 
        continue, 
    end
    
    ext = exts(it).extent;
	for ipoint = ipoints(i_from_ends2center)  %#ok<ALIGN>
        
        for xdir = [-1 1] %#ok<ALIGN>
%            if xdir==-1, exts(it).halign='Right'; else exts(it).halign='Left'; end
            exts(it).halign='Left';
            
            for ydir = [-1 1] %#ok<ALIGN>
%                if ydir==-1, exts(it).valign='Top'; else exts(it).valign='Bottom'; end
                exts(it).valign='Bottom';

                xbox(1) = x(ipoint);
                xbox(2) = xbox(1) + ext(3)*xdir;
                ybox(1) = y(it,ipoint) + ymgn_ownline*ydir; %assuming exts and y have the same length
                ybox(2) = ybox(1) + ext(4)*ydir;
                xbox = sort(xbox);
                ybox = sort(ybox);

                exts(it).xstart = xbox(1);
                exts(it).xend = xbox(2);
                exts(it).ystart = ybox(1);
                exts(it).yend = ybox(2);

                if isnan(exts(it).ystart) || isnan(exts(it).yend)
                    continue;
                end
                
                if exts(it).xstart<xlimits(1)+xmgn_box || exts(it).xend>xlimits(2)-xmgn_box || ... %check crossing picture borders
                   exts(it).ystart<ylimits(1)+ymgn_box || exts(it).yend>ylimits(2)-ymgn_box 
                    continue;
                end

                if fl_debug
                     garbobjs(end+1) = text(exts(it).xstart, exts(it).ystart, exts(it).str, 'color','r', 'FontSize', exts(it).fontsz, ...
                         'HorizontalAlignment',exts(it).halign, 'VerticalAlignment',exts(it).valign, ...
                         'EdgeColor','r','Margin',eps);
                     plot(exts(it).xstart, exts(it).ystart, 'dr', 'MarkerFaceColor','r' );
                end

                ok = 1;
                for ip=1:numel(buf.h_plots)
                    if exts(ip).disabled 
                        continue
                    end
                    
                    if ip==it %for own line do not use border
                        ix_in_ext = find(x>=exts(it).xstart & x<=exts(it).xend);
                        if any(y(ip,ix_in_ext)>=exts(it).ystart & y(ip,ix_in_ext)<=exts(it).yend)
                            ok = 0; %%extent overlapes line
                            break;
                        end
                    else %using borders for alien lines
                        ix_in_ext = find(x>=exts(it).xstart-xmgn_alien & x<=exts(it).xend+xmgn_alien);
                        if any(y(ip,ix_in_ext)>=exts(it).ystart-ymgn_alien & y(ip,ix_in_ext)<=exts(it).yend+ymgn_alien)
                            ok = 0; %%extent overlapes line
                            break;
                        end
                    end
                end
                if ~ok
                    continue;            
                end

                ok = 1;
                for it2=1:it-1 %check overlaping with prev organized extents
                    if exts(it2).disabled, continue, end;
                    if rectint([exts(it2).xstart exts(it2).ystart exts(it2).extent(3) exts(it2).extent(4)],...
                               [exts(it).xstart exts(it).ystart exts(it).extent(3) exts(it).extent(4)] ) > 0 %check intersection area
                        ok = 0;
                        break;
                    end
                end
                if ok
                    exts(it).success = 1;
                    if fl_debug
                        set(garbobjs(end),'Color','g');
                    end
                    res = res + 1;
                    break;
                end

            end %ydir = [-1 1]
            if exts(it).success, break, end
            
        end %xdir = [-1 1]
        if exts(it).success, break, end
        
    end %ipoint = ipoints 
end %it=1:numel(exts)
disp(['Successfully moved:' int2str([exts.success])])

for it=1:numel(buf.ht)
    if (exts(it).disabled == 1) 
        if ~flags.develmode %in devel mode move nothing
            disp(['Label for similar plot deleted: ' exts(it).str])
            delete(exts(it).ht);
            disp(['similar plot deleted: ' exts(it).str])
            delete(exts(it).h_plots);
        end
    else
        set(exts(it).ht,'Position',[exts(it).xstart exts(it).ystart 0], ...
            'HorizontalAlignment',exts(it).halign, 'VerticalAlignment',exts(it).valign);
    end
end
if fl_debug
    delete(garbobjs)
end
