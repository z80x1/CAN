%convert r9 (Cyd) unique conformations to  2'deoxy-5methylcytidine (r270) ones by 
%replacing H5 atom with C52 and addind H51, H52 and H53 at 1A
%
% Version 1.0    
% Last modified  R O Zhurakivsky 2008-05-19
% Created        R O Zhurakivsky 2006-?-?

clear 
format compact
atomsind

%--------------------------------------
moltype = 9  %#ok
usedpackage='Gaussian'  %#ok
theory='dftV2'  %#ok
onlyoriginal=1  %#ok %process db with only original conformations
workname='r27001' %#ok
%--------------------------------------


workdbname=['r' int2str(moltype)]   %#ok
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


    %changing bH5 to bCC5 and placing at needed distance
    indA=find(ms0.pind==find(strcmp(pind.labels,'bH5')));
    ms0.pind(indA) = find(strcmp(pind.labels,'bCC5'));
    ms0.labels(indA)={'C'};

    distCC=1.5; %C5-CC5 bond distance from optimized thymidine structure
    indB=find(ms0.pind==find(strcmp(pind.labels,'bC5')));
    A=[ms0.x(indA) ms0.y(indA) ms0.z(indA)];
    B=[ms0.x(indB) ms0.y(indB) ms0.z(indB)];
    dirC5C52=A-B;
    A=B+dirC5C52/norm(dirC5C52)*distCC;
    ms0.x(indA)=A(1);
    ms0.y(indA)=A(2);
    ms0.z(indA)=A(3);


    %adding bH51, bH52 and bH53
    indA=find(ms0.pind==find(strcmp(pind.labels,'bC4')));
    indB=find(ms0.pind==find(strcmp(pind.labels,'bC5')));
    indC=find(ms0.pind==find(strcmp(pind.labels,'bC6')));
    indD=find(ms0.pind==find(strcmp(pind.labels,'bCC5')));
    A=[ms0.x(indA) ms0.y(indA) ms0.z(indA)];
    B=[ms0.x(indB) ms0.y(indB) ms0.z(indB)];
    C=[ms0.x(indC) ms0.y(indC) ms0.z(indC)];
    D=[ms0.x(indD) ms0.y(indD) ms0.z(indD)];

    %place H51 and H52 in C4C5C6 plane and so that C52H51 is parallel to C6H5 and C52H52 is parallel to C4C5
    % C52H53 is parpendicular both to C52H51 and C52H52 one	(CH distance is equal to 1.1A)
    distCH=1.095;
    andCCO=111.36;

    C6C5=B-C;
	C4C5=B-A;
	xH51=D+C6C5/norm(C6C5)*distCH;
	xH52=D+C4C5/norm(C4C5)*distCH;
	dirC52H53=cross(C6C5,C4C5);
	xH53=D+dirC52H53/norm(dirC52H53)*distCH;

	C5C52=D-B;


    ms0.atomnum=ms0.atomnum+3;
    ms0.x=[ms0.x; xH51(1); xH52(1); xH53(1)];
    ms0.y=[ms0.y; xH51(2); xH52(2); xH53(2)];
    ms0.z=[ms0.z; xH51(3); xH52(3); xH53(3)];
    ms0.labels=[ms0.labels; {'H'}; {'H'}; {'H'}];
    ms0.ind=[ms0.ind; max(ms0.ind)+1; max(ms0.ind)+2; max(ms0.ind)+3];
    ms0.pind(end+1)=find(strcmp(pind.labels,'bH51'));
    ms0.pind(end+1)=find(strcmp(pind.labels,'bH52'));
    ms0.pind(end+1)=find(strcmp(pind.labels,'bH53'));


    ms0=createbondtable(ms0);

    ms0.desc=[workname,'_',sprintf('%03d',fileind),'_',ms0.desc];

    order=[6 1 2]; 

    or = [order, setdiff(1:ms0.atomnum,order)]; %!!if using see rot14 for better inplementation
    savemol(odir,ms0,0,or);

%not needed
%    [xxx,order]=sortrows(ms0.pind);

%   savemolgs(odir,ms0,3,order,gtemplname); %Gaussian with XYZ
    savemolgs(odir,ms0,4,order,fullgtemplname); %Gaussian with ZMT

end


