%convert AC6 O2 protonated unique conformations by adding H3 at 1.084A from N3
%
% Version 1.0    
% Last modified  R O Zhurakivsky 2008-05-19
% Created        R O Zhurakivsky 2012-07-30

clear 
format compact
atomsind
pindsdef

workdbname=[CD.dbdir filesep 'r14_g_dftV3bis_or.mat']   %#ok
workname='r24101'   %#ok
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

    indA=find(ms0.pind==find(strcmp(pind.labels,'bC2')));
    indB=find(ms0.pind==find(strcmp(pind.labels,'bN3')));
    indC=find(ms0.pind==find(strcmp(pind.labels,'bC4')));
    A=[ms0.x(indA) ms0.y(indA) ms0.z(indA)];
    B=[ms0.x(indB) ms0.y(indB) ms0.z(indB)];
    C=[ms0.x(indC) ms0.y(indC) ms0.z(indC)];
    X=(A+C)./2;

    ortXB=(B-X)/norm(B-X);
    
    dist=1.084;
    D=B+dist*ortXB;

    ms0.atomnum=ms0.atomnum+1;
    ms0.x(end+1)=D(1);
    ms0.y(end+1)=D(2);
    ms0.z(end+1)=D(3);
    ms0.labels(end+1)={'H'};
    ms0.ind(end+1)=max(ms0.ind)+1;
    ms0.pind(end+1)=find(strcmp(pind.labels,'bH3'));
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