%Isakova Alina
%r15 (dAdo) to r150 (dAdo imino)
%
% Version 1.0    
% Last modified  R O Zhurakivsky 2008-05-19
% Created        R O Zhurakivsky 2006-?-?


clear 
format compact
atomsind

workdbname=[CD.dbdir filesep 'r15_g_dftV2_or.mat']    %#ok
workname='r15001' %#ok
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

    %changing bH2 to bCN2. C-N distance equals 1.35
    inde=find(ms0.pind==find(strcmp(pind.labels,'bH2')));
    ms0.pind(inde) = find(strcmp(pind.labels,'bCN2'));
    ms0.labels(inde)={'N'};
    
    distCN=1.35; %CN distance is equal 1.35
    
    indf=find(ms0.pind==find(strcmp(pind.labels,'bC2')));
    
    e=[ms0.x(inde) ms0.y(inde) ms0.z(inde)];
    f=[ms0.x(indf) ms0.y(indf) ms0.z(indf)];

    ef = e-f;
    e=f+ef/norm(ef)*distCN;
    
    ms0.x(inde)=e(1);
    ms0.y(inde)=e(2);
    ms0.z(inde)=e(3);
    
    %adding H21 and H22
    indA=find(ms0.pind==find(strcmp(pind.labels,'bN1')));
    indB=find(ms0.pind==find(strcmp(pind.labels,'bC2')));
    indC=find(ms0.pind==find(strcmp(pind.labels,'bN3')));
    indD=find(ms0.pind==find(strcmp(pind.labels,'bCN2')));
    A=[ms0.x(indA) ms0.y(indA) ms0.z(indA)];
    B=[ms0.x(indB) ms0.y(indB) ms0.z(indB)];
    C=[ms0.x(indC) ms0.y(indC) ms0.z(indC)];
    D=[ms0.x(indD) ms0.y(indD) ms0.z(indD)];
    
    %place H21 and H22 in N1C2N3 plane so that CN2H21 is parallel to N3C2 and CN2H22 is parallel to N1C2
    distNH=1.01;
    
    C2N3=B-C;
	N1C2=B-A;
	xH21=D+C2N3/norm(C2N3)*distNH;
	xH22=D+N1C2/norm(N1C2)*distNH;
	
    ms0.atomnum=ms0.atomnum+2;
    ms0.x=[ms0.x; xH21(1); xH22(1)];
    ms0.y=[ms0.y; xH21(2); xH22(2)];
    ms0.z=[ms0.z; xH21(3); xH22(3)];
    ms0.labels=[ms0.labels; {'H'}; {'H'}];
    ms0.ind(end+1)= max(ms0.ind)+1; 
    ms0.ind(end+1)= max(ms0.ind)+1;
    ms0.pind(end+1)=find(strcmp(pind.labels,'bH21'));
    ms0.pind(end+1)=find(strcmp(pind.labels,'bH22'));
 
    
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
    
   
    ms0=createbondtable(ms0);

    ms0.desc=[workname,'_',sprintf('%03d',fileind),'_',ms0.desc];

    order=[6 1 2]; 

    or = [order, setdiff(1:ms0.atomnum,order)]; %!!if using see rot14 for better inplementation
    savemol(odir,ms0,0,or);

%not needed
%    [xxx,order]=sortrows(ms0.pind);
   savemolgs(odir,ms0,3,order,fullgtemplname); %Gaussian with XYZ
%    savemolgs(odir,ms0,4,order,fullgtemplname); %Gaussian with ZMT

end


