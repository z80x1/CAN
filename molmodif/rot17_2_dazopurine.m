%convert r15 (dAdo) unique conformations to  deoxy8azopurine (r162) ones by 
%1)deleting NH2 group 
%2) changing C8H to N8
%
% Version 1.0    
% Last modified  R O Zhurakivsky 2008-05-19
% Created        R O Zhurakivsky 2006-?-?

format compact
atomsind

workdbname=[CD.dbdir filesep 'r15_g_dftV2_or.mat']    %#ok
workname='r16201' %#ok
sortmode=2      %#ok

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

    
     %changing bCN6 to bH6
    indp=find(ms0.pind==find(strcmp(pind.labels,'bCN6')));
    ms0.pind(indp) = find(strcmp(pind.labels,'bH6'));
    ms0.labels(indp)={'H'};

    distCH=1.095; %C-H distance is equal 1.095
    
    indq=find(ms0.pind==find(strcmp(pind.labels,'bC6')));
    
    p=[ms0.x(indp) ms0.y(indp) ms0.z(indp)];
    q=[ms0.x(indq) ms0.y(indq) ms0.z(indq)];

    pq = p-q;
    p=q+pq/norm(pq)*distCH;
    
    ms0.x(indp)=p(1);
    ms0.y(indp)=p(2);
    ms0.z(indp)=p(3);
        
    %deleting H61 and H62
    ind=find(ms0.pind==find(strcmp(pind.labels,'bH61')));
    ms0=delatom2(ms0,ind);
    ind=find(ms0.pind==find(strcmp(pind.labels,'bH62')));
    ms0=delatom2(ms0,ind);

    %changing bC8 to bN8
    indp=find(ms0.pind==find(strcmp(pind.labels,'bC8')));
    ms0.pind(indp) = find(strcmp(pind.labels,'bN8'));
    ms0.labels(indp)={'N'};
    %deleting H8
    ind=find(ms0.pind==find(strcmp(pind.labels,'bH8')));
    ms0=delatom2(ms0,ind);
    
   
    ms0=createbondtable(ms0);

    ms0.desc=[workname,'_',sprintf('%03d',fileind),'_',ms0.desc];

    order=[6 1 2]; 

    or = [order, setdiff(1:ms0.atomnum,order)]; %!!if using see rot14 for better implementation
    savemol(odir,ms0,0,or);

%not needed
%    [xxx,order]=sortrows(ms0.pind);
   savemolgs(odir,ms0,3,order,fullgtemplname); %Gaussian with XYZ
%    savemolgs(odir,ms0,4,order,fullgtemplname); %Gaussian with ZMT

end


