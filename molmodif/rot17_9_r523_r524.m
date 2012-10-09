%convert d5ClU unique conformations to d5FU ones by changing 5Cl to 5F
%
% Version 1.0    
% Created        R O Zhurakivsky 2010-02-17

clear 
format compact
atomsind
pindsdef

workdbname=[CD.dbdir filesep 'r523_g_dftV3.mat']   %#ok
workname='r52401'   %#ok
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

    ibCl5 = find(ms0.pind==find(strcmp(pind.labels,'bCl5')));
    ibC5 = find(ms0.pind==find(strcmp(pind.labels,'bC5')));
    ms0.labels(ibCl5)={'F'};

    A=[ms0.x(ibC5) ms0.y(ibC5) ms0.z(ibC5)];
    B=[ms0.x(ibCl5) ms0.y(ibCl5) ms0.z(ibCl5)];
    ortAB=(B-A)/norm(B-A);
    dist=1.4;
    C=A+dist*ortAB;
    ms0.x(ibCl5)=C(1);
    ms0.y(ibCl5)=C(2);
    ms0.z(ibCl5)=C(3);

    ms0.pind(ibCl5)=find(strcmp(pind.labels,'bF5'));

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


