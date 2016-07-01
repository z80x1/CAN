%rot37_Hbonddynamic:    Plot energy of H-bond while molecule is rotating e.g. round chi
%
% Version 1.0    
% Last modified  R O Zhurakivsky 2009-10-20
% Created        R O Zhurakivsky 2009-05-25

tic
clear 
format compact

atomsind
pindsdef

%---------------------------------
moltype=[9 13 15 16] %#ok
usedpackage='Gaussian' %#ok
theory='dftV3' %#ok
onlyoriginal=0;  % process db with only original conformations
molconf='EabcA';
onlytruehbonds = 2 %#ok 1 - only 3 or 4 atoms Hbonds, 2 - only chosen Hbonds

chosenhbonds = [{'bC6bH6...pO5'} {'bC8bH8...pO5'} {'pC1pH12...bO2'} {'bC6bH6...pO4'} {'pC2pH21...pO5'}];

%bondstr='bC6bH6...pO4';
%---------------------------------
pcolor=[{'r'} {'g'} {'b'} {'k'} {'m'} {'c'} ...
        {[.8627 .0784 .2353]} {[.5412 .1686 .8863]}  ...
        {[0 .5 0]} {[.5 0 0]} {[0 0 .5]} {[.5 .5 0]} {[.5 0 .5]} {[0 .5 .5]} {[.5 .5 .5]} ...
        {[.5 .1 .9]} {[.8 .2 .2]} {[.8 .8 .2]}...
        {[.9 .4 .9]} {[.2 .4 .6]} {[.6 .4 .6]} {[.6 .2 .2]} {[.8 .2 .8]} ...
        {[.2 .8 .8]} {[.2 .8 .2]} {[.2 .2 .8]} {[.4 .9 .1]} {[.1 .3 .6]} {'y'} {[.75 .75 .75]} {[.2745 .5098 .7059]}];
%{[1 .2706 0]}  - almost red
%psign='x+*osdv^<>ph';
psign='*d+xv^<>hp.o';

allbondsstr = [];
for mind=1:numel(moltype)

    workdbname0=['r' int2str(moltype(mind))] %#ok
    if ~isempty(molconf)
        workdbname0=[workdbname0 '_' molconf];
    end
    workdbname=workdbname0;
    if strcmp(usedpackage,'Gaussian')
      workdbname=[workdbname '_g'];
    end
    if ~strcmp(theory,'dft')
      workdbname=[workdbname '_' theory];
    end
    if onlyoriginal
        templ='_or';
        workdbname = [workdbname templ];
    end
    workdbname=[CD.dbdir filesep workdbname '.mat'] %#ok



    if exist(workdbname,'file')
        load(workdbname,'workdb')
    else
        workdb=[];
    end
    recnum=numel(workdb);

    if isfield(workdb,'AIM')
        AIM=[workdb.AIM];
    else
        AIM.desc=[];
    end

     
    if onlytruehbonds ~= 2
        if onlytruehbonds == 0
            bondsstr=unique([AIM.desc]);
        elseif onlytruehbonds == 1
            bondsstr0=unique([AIM.desc]);
            bondsstr={};
            for ii=1:numel(bondsstr0)
                if numel(bondsstr0{ii})>=12 %only bonds with 3 atoms or more
                    bondsstr(end+1)=bondsstr0(ii);
                end
            end
        else
            error('incorrect onlytruehbonds value');
        end
        allbondsstr = [allbondsstr setdiff(bondsstr,allbondsstr)];
    else %only chosen Hbonds
        bondsstr = chosenhbonds;
        allbondsstr = chosenhbonds;
    end
        

    for bind=1:numel(bondsstr)
        
        energy=[];
        param=[];
        for i=1:recnum

            sdesc{i} = workdb(i).prop.sdesc;
            if isempty(workdb(i).AIM), continue, end;
            ind=strcmpcellar(workdb(i).AIM.desc,bondsstr(bind));
            if ind
               param(end+1) = workdb(i).param;
               energy(end+1) = abs(0.5*workdb(i).AIM.V(ind)*CC.encoef);

            end
        end
        
        allbind = strcmpcellar(allbondsstr,bondsstr(bind));
        allparam(mind,allbind)={param};
        allenergy(mind,allbind)={energy};
    end
end %mind


    h=[]; %plot handles array
    legs=[]; %legends array
    fl_firstplot=1;
    for mind=1:size(allparam,1)
    for bind=1:size(allparam,2)
        if fl_firstplot
            figure
        end
        if ~isempty(allparam{mind,bind})
            if mind==4 %crosses are too small - lets increase their size 
                markersize=6;
            else
                markersize=4;
            end
            h(end+1)=plot(allparam{mind,bind},allenergy{mind,bind},psign(mod(mind-1,numel(psign))+1),'Color', pcolor{bind},'MarkerSize', markersize);
            legs{end+1}=['r' int2str(moltype(mind)) ' ' allbondsstr{bind} ];
        end
        if fl_firstplot
            hold on
            fl_firstplot=0;
        end
    end
    end

    hl=legend(h, legs,'Location','NEO');
    set(hl,'FontSize',6 );
    add = '';
    if onlytruehbonds==1
        add = ' (w/o 2 atoms contacts)';
    elseif onlytruehbonds==2
        add = ' (only chosen)';
    end
        
    title(['r' int2str(moltype) ' ' molconf ' H-bonds energy of chi dependence' add])
    grid on
    xlabel('\chi');
    ylabel('bond energy, kcal/mol');

    oldaxis=axis;
    axis([-180 180 0 oldaxis(4)]);
    hold off

if 0
    param=[];
    GOenergy=[];
    Henergyall=[];
    for i=1:recnum
        if ~isempty(workdb(i).AIM)
            Henergy=0;
            for hind=1:numel(workdb(i).AIM.desc)
    if 0
                if strcmp(workdb(i).AIM.desc(hind),'pO4...bO2') 
                    
                elseif strcmp(workdb(i).AIM.desc(hind),'pO5...bN1')
                    Henergy=Henergy-abs(0.5*workdb(i).AIM.V(hind)*CC.encoef);
                elseif strcmp(workdb(i).AIM.desc(hind),'pO5...bO2') 
                    Henergy=Henergy-abs(0.5*workdb(i).AIM.V(hind)*CC.encoef);
                elseif strcmp(workdb(i).AIM.desc(hind),'pC2...bO2') 
                    
                elseif strcmp(workdb(i).AIM.desc(hind),'pC2pH21...pO5')
                    
                elseif strcmp(workdb(i).AIM.desc(hind),'pC2pH21...bH6bC6')
                    
                else
                    Henergy=Henergy+abs(0.5*workdb(i).AIM.V(hind)*CC.encoef);
                end
    else
                    Henergy=Henergy+abs(0.5*workdb(i).AIM.V(hind)*CC.encoef);
    end
            end
            Henergyall(end+1)=Henergy;
        else
            continue
        end
        param(end+1)=workdb(i).param;
        GOenergy(end+1)=workdb(i).GO.energy;
    end

    if 0
        figure
        hold on
        GOenergy=(GOenergy-min(GOenergy))*CC.encoef;
        Henergyall=max(Henergyall)-Henergyall;
        [param,I]=sort(param);
        GOenergy=GOenergy(I);
        Henergyall=Henergyall(I);
        h1=plot(param,GOenergy,'r.-');
        h2=plot(param,Henergyall,'b.-');
        title([strrep(workdbname0,'_',' ') ' Sum H-bonds energy']);
        xlabel('\chi, degree');
        ylabel('bond energy, kcal/mol');
        oldaxis=axis;
        axis([-180 180 0 oldaxis(4)]);
        grid on
        hl=legend([h1; h2], [{'GOenergy'}; {'Henergyall'}],'Location','NEO');
        hold off
    end
end
toc

if 1 %Hbonds existance ranges
    for mind=1:size(allparam,1)
    for bind=1:size(allparam,2)
        a=allparam{mind,bind};
        if ~isempty(a)
            a1 = a;
            a1(find(a1<0))=a1(find(a1<0))+360;
            if (max(a)-min(a))+1 > (max(a1)-min(a1))
                a=a1;
            end
            ['r' int2str(moltype(mind)) ' ' allbondsstr{bind} ',min' int2str(min(a))  ',max' int2str(max(a)) ]
        end
    end
    end
end