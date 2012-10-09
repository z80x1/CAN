function pl2ax(plSource,axDest,addtitle,ky)
% Copy plot to another axes.
% Still doesn't do legend. Don't know abouut colorbars, but hey its a start.
if nargin<4
    ky=1;
end
if nargin<3
    addtitle='';
elseif ~isempty(addtitle)
    addtitle=[': ' addtitle];
end


axSource = get(plSource,'Parent');
% Find children, copy non-text ones.
kids=allchild(axSource);
nontextkids=kids(~strcmp(get(kids,'type'),'text'));
%copyobj(nontextkids,axDest);
plDest=copyobj(plSource,axDest);
set(plDest,'Marker','.','MarkerSize',6);

% Axes Directions (may need to add other properties)
meth={'YDir','XLim','XGrid','XMinorGrid','XTick','XTickLabel'};
cellfun(@(m)set(axDest,m,get(axSource,m)),meth);

% Special treatment for text-children
%subplot(hDest)
xlabel(axDest,get(get(axSource,'xlabel'),'string'));
ylabel(axDest,get(get(axSource,'ylabel'),'string'));
title(axDest,[get(get(axSource,'title'),'string') addtitle],'FontSize',9);

old_axis=axis(axDest);
y=get(plDest,'YData');
ysize=max(y)-min(y);
axis(axDest,[old_axis(1) old_axis(2) min(y)-0.01*ysize max(y)+0.01*ysize]);

set(axDest,'XGrid','on','YGrid','on','GridLineStyle','--');
%set(axDest,'XMinorGrid','on','YMinorGrid','on','MinorGridLineStyle',':');


end 