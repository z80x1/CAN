%change r161 (deoxypurine) ro r15 (dAdo) 
%by changing H6 to bCN6 and adding bH61,bH62 atoms 
%
% Version 1.0    
% Last modified  R O Zhurakivsky 2008-05-19
% Created        R O Zhurakivsky 2006-?-?

clear 
format compact
atomsind
pindsdef

%-------------------------------------
workdbname=[CD.dbdir filesep 'r161_g_dftV3_or.mat']    %#ok
workname='rf30' %#ok
sortmode=2      %#ok

odir=[CD.xyzdir filesep workname];

includelist.do=1;   %#ok
includelist.type='file';
includelist.name=[odir filesep 'includelist'];
includelist.rectype = 'sdesc' %may have 'sdesc' or 'index' values
includelist  %#ok
%-------------------------------------


global pind
load(workdbname,'workdb')

gtemplname=[workname '_templ.gjf']  %#ok
fullgtemplname=[CD.templatesdir filesep gtemplname];
odir=[CD.xyzdir filesep workname];
if exist(odir,'dir')~=7
   mkdir(odir);
end

inclset={};
inclsetind=[];
if includelist.do
if ~isempty(includelist.name)

  if strcmp(includelist.type,'dir') %choose record with identifier in pointed directory
    inclfiles = dir(includelist.name);
    numfiles = size(inclfiles,1);
  elseif strcmp(includelist.type,'file')

    [fid,emessage]=fopen(includelist.name);
    if fid==-1
      error(['No file ' includelist.name ' with descriptions to include found: ' emessage]);
    end

    tline=fgetl(fid);
    i=1;
    while tline~=-1
      inclfiles(i).name=tline;
      tline=fgetl(fid);
      i=i+1;
    end 
    fclose(fid);
    numfiles=i-1;

  end

  if ~numfiles
    disp('No descriptions found to include.')
  else
    for i=1:numfiles
      if strcmp(includelist.rectype,'sdesc')

          dlm1=strfind(inclfiles(i).name,'_');
          if isempty(dlm1) 
            dlm1=0;
          end   
          dlm2=strfind(inclfiles(i).name,'.');
          if isempty(dlm2), dlm2=size(inclfiles(i).name,2)+1; end;
          inclset(i) = {inclfiles(i).name((dlm1(end)+1):(dlm2-1))}; %complete set of names
      else
          inclsetind(i) = sscanf(inclfiles(i).name,'%d'); %complete set of indexes of records to include
      end
    end
  end
end
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
    if ~isempty(inclset) && isempty(strcmpcellar(inclset,workdb(i).prop.sdesc)), continue, end;

    fileind=fileind+1;

    ms0.labels=workdb(i).labels;
    ms0.x=workdb(i).x;
    ms0.y=workdb(i).y;
    ms0.z=workdb(i).z;
    ms0.atomnum=workdb(i).atomnum;
    ms0.ind=workdb(i).ind;
    ms0.desc=workdb(i).prop.sdesc;
    ms0.pind=workdb(i).pind;

    %changing bH6 to bCN6. C-N distance equals 1.35
    inde=find(ms0.pind==find(strcmp(pind.labels,'bH6')));
    ms0.pind(inde) = find(strcmp(pind.labels,'bCN6'));
    ms0.labels(inde)={'N'};
    
    distCN=1.35; %CN distance is equal 1.35
    
    indf=find(ms0.pind==find(strcmp(pind.labels,'bC6')));
    
    e=[ms0.x(inde) ms0.y(inde) ms0.z(inde)];
    f=[ms0.x(indf) ms0.y(indf) ms0.z(indf)];

    ef = e-f;
    e=f+ef/norm(ef)*distCN;
    
    ms0.x(inde)=e(1);
    ms0.y(inde)=e(2);
    ms0.z(inde)=e(3);
    
    %adding H61 and H62
    indA=find(ms0.pind==find(strcmp(pind.labels,'bN1')));
    indB=find(ms0.pind==find(strcmp(pind.labels,'bC6')));
    indC=find(ms0.pind==find(strcmp(pind.labels,'bC5')));
    indD=find(ms0.pind==find(strcmp(pind.labels,'bCN6')));
    A=[ms0.x(indA) ms0.y(indA) ms0.z(indA)];
    B=[ms0.x(indB) ms0.y(indB) ms0.z(indB)];
    C=[ms0.x(indC) ms0.y(indC) ms0.z(indC)];
    D=[ms0.x(indD) ms0.y(indD) ms0.z(indD)];
    
    %place H61 and H62 in N1C6C5 plane so that CN6H61 is parallel to C5C6 and CN6H62 is parallel to N1C6
    distNH=1.01;
    
    BC=B-C;
	BA=B-A;
	xH1=D+BC/norm(BC)*distNH;
	xH2=D+BA/norm(BA)*distNH;
	
    ms0.atomnum=ms0.atomnum+2;
    ms0.x=[ms0.x; xH1(1); xH2(1)];
    ms0.y=[ms0.y; xH1(2); xH2(2)];
    ms0.z=[ms0.z; xH1(3); xH2(3)];
    ms0.labels=[ms0.labels, {'H'}, {'H'}];
    ms0.ind(end+1)= max(ms0.ind)+1; 
    ms0.ind(end+1)= max(ms0.ind)+1;
    ms0.pind(end+1)=find(strcmp(pind.labels,'bH61'));
    ms0.pind(end+1)=find(strcmp(pind.labels,'bH62'));
 
    
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


