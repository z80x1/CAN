%convert d5BrU unique conformations to d5ClU ones by changing 5Br to 5Cl
%
% Version 1.0    
% Created        R O Zhurakivsky 2009-11-19

clear 
format compact
atomsind
pindsdef

workdbname=[CD.dbdir filesep 'r521_g_dftV3.mat']   %#ok
workname='r52301'   %#ok
sortmode=2          %#ok

global pind
load(workdbname,'workdb')

gtemplname=[workname '_templ.gjf']  %#ok
fullgtemplname=[CD.templatesdir filesep gtemplname];
odir=[CD.xyzdir filesep workname];
if exist(odir,'dir')~=7
   mkdir(odir);
end

GOenergy=[];
sdesc={};
recnum=numel(workdb);
for i=1:recnum

    if isfield(workdb(i),'GO') && isfield(workdb(i).GO,'energy')
        GOenergy(i) = workdb(i).GO.energy;
    else
        sortbyenergy=0;
    end
    sdesc(i) = {workdb(i).prop.sdesc};
end

if sortmode==1 %sort by sdesc
  [xxx,sortind]=sort(sdesc);
elseif sortmode==2 %sort by energy
  [xxx,sortind]=sort(GOenergy);
else %without sorting
  sortind=1:recnum;
end

fileind=0;
for i=sortind

    if workdb(i).new~='Y'
       continue
    end  
    fileind=fileind+1;

    ms0.labels=workdb(i).labels;
    ms0.x=workdb(i).x;
    ms0.y=workdb(i).y;
    ms0.z=workdb(i).z;
    ms0.atomnum=workdb(i).atomnum;
    ms0.ind=workdb(i).ind;
    ms0.desc=workdb(i).prop.sdesc;
    ms0.pind=workdb(i).pind;

    ind=find(ms0.pind==find(strcmp(pind.labels,'bBr5')));
    ms0.labels(ind)={'Cl'};

    ms0=createbondtable(ms0);

    ms0.desc=[workname,'_',sprintf('%03d',fileind),'_',ms0.desc];

    pinds=[6 1 2]; %pinds pO4, pC1, pC2
    order=[];
    for I=1:numel(pinds)
      order(end+1)=find(ms0.pind==pinds(I));
    end
    [xxx,or]=setdiff(ms0.pind,pinds); %output xyz with chosen atoms first then in common way - sorted by pind

    savemol(odir,ms0,0,or);
    savemolgs(odir,ms0,4,order,fullgtemplname); %Gaussian with ZMT

end


