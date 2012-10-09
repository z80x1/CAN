%normal pirimidins (Cyd,Urd,Thd) to 5AZA or 6AZA conversion
%e.g.: convert r9 (Cyd) unique conformations to  5azaCyd (rl/r21) by deleting bH5 atom and changing bC5 to bN5
%
% Version 1.0    
% Last modified  R O Zhurakivsky 2008-05-19
% Created        R O Zhurakivsky 2006-?-?

clear 
format compact
atomsind

%workdbname=[CD.dbdir filesep 'r9_g_dftV2_or.mat']
workdbname=[CD.dbdir filesep 'r12_g_dftV2_or.mat']
%workname='rl04'
%workname='rh04'
workname='rj04'
sortmode=2

fl_to5aza=1
fl_to6aza=0

global pind
load(workdbname,'workdb')

gtemplname=[workname '_templ.gjf']
fullgtemplname=[CD.templatesdir filesep gtemplname];
odir=[CD.xyzdir filesep workname];
if exist(odir)~=7
   mkdir(odir);
end

recnum=numel(workdb);
for i=1:recnum

    if isfield(workdb(i),'GO') & isfield(workdb(i).GO,'energy')
        GOenergy(i) = workdb(i).GO.energy;
    else
        sortbyenergy=0;
    end
    sdesc(i) = {workdb(i).prop.sdesc};
end

if sortmode==1 %sort by sdesc
  [xxx,ind]=sort(sdesc);
elseif sortmode==2 %sort by energy
  [xxx,ind]=sort(GOenergy);
else %without sorting
  ind=1:recnum;
end

fileind=0;
for i=ind

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

    if fl_to5aza
       %changing bC5 to bN5
       indA=find(ms0.pind==find(strcmp(pind.labels,'bC5')));
       ms0.pind(indA) = find(strcmp(pind.labels,'bN5'));
       ms0.labels(indA)={'N'};

       %deleting bH5
       ind=find(ms0.pind==find(strcmp(pind.labels,'bH5')));
       ms0=delatom2(ms0,ind);
    end
    if fl_to6aza
       %changing bC6 to bN6
       indA=find(ms0.pind==find(strcmp(pind.labels,'bC6')));
       ms0.pind(indA) = find(strcmp(pind.labels,'bN6'));
       ms0.labels(indA)={'N'};

       %deleting bH6
       ind=find(ms0.pind==find(strcmp(pind.labels,'bH6')));
       ms0=delatom2(ms0,ind);
    end

    ms0=createbondtable(ms0);

    ms0.desc=[workname,'_',sprintf('%03d',fileind),'_',ms0.desc];

    pinds=[6 1 2]; %pinds pO4, pC1, pC2
    order=[];
    for I=1:numel(pinds)
      order(end+1)=find(ms0.pind==pinds(I));
    end
    [xxx,or]=setdiff(ms0.pind,pinds); %output xyz with chosen atoms first then in common way - sorted by pind
%not needed
%    %ATTENTION: setdiff outputs row vector with sorted values
%    order=[order order1];

    savemol(odir,ms0,0,or);
    savemolgs(odir,ms0,4,order,fullgtemplname); %Gaussian with ZMT

end


