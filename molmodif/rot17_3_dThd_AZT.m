%convert r13 (dThd) unique conformations to  AZT (r413) ones by replacing O3'H group with N-N-N group 
%
% Version 1.0    
% Last modified  R O Zhurakivsky 2010-10-15
% Created        R O Zhurakivsky 2010-10-15

clear 
format compact
atomsind
pindsdef

workdbname=[CD.dbdir filesep 'r13_g_dftV3_or.mat']    %#ok
workname='r41301' %#ok
sortmode=2      %#ok
fl_save = 1;

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



    
    %adding N-N-N
    indA=find(ms0.pind==find(strcmp(pind.labels,'pO3')));
    ms0.pind(indA) = find(strcmp(pind.labels,'pN31'));
    ms0.labels(indA)={'N'};

    indB=find(ms0.pind==find(strcmp(pind.labels,'pH32')));
    ms0.pind(indB) = find(strcmp(pind.labels,'pN32'));
    ms0.labels(indB)={'N'};

    A=[ms0.x(indA) ms0.y(indA) ms0.z(indA)];
    B=[ms0.x(indB) ms0.y(indB) ms0.z(indB)];
    ortAB=(B-A)/norm(B-A);
    distNN=1.1;
    Bbis=A+distNN*ortAB;
    ms0.x(indB)=Bbis(1);
    ms0.y(indB)=Bbis(2);
    ms0.z(indB)=Bbis(3);


    C=Bbis+distNN*ortAB;

    ms0.atomnum=ms0.atomnum+1;
    ms0.x(end+1)=C(1);
    ms0.y(end+1)=C(2);
    ms0.z(end+1)=C(3);
    ms0.labels(end+1)={'N'};
    ms0.ind(end+1)=max(ms0.ind)+1;
    ms0.pind(end+1)=find(strcmp(pind.labels,'pN33'));


    ms0=createbondtable(ms0);

    ms0.desc=[workname,'_',sprintf('%03d',fileind),'_',ms0.desc];

    order=[6 1 2]; 

    or = [order, setdiff(1:ms0.atomnum,order)]; %!!if using see rot14 for better inplementation

%plotmol(ms0,'r',0);
%keyboard
    if fl_save
      savemol(odir,ms0,0,or);

%not needed
%    [xxx,order]=sortrows(ms0.pind);
%     savemolgs(odir,ms0,3,order,gtemplname); %Gaussian with XYZ
      savemolgs(odir,ms0,4,order,fullgtemplname); %Gaussian with ZMT
    end

end


