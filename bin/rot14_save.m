%rot14: save geomerties and/or create input files from database data
%
% Version 1.0    
% Last modified  R O Zhurakivsky 2008-08-01
% Created        R O Zhurakivsky 2005-09-?

%051122 aDded includelist processing

format compact
clear 
%clear sdesc GOenergy
atomsind

%--------------------------
moltype=9  %#ok
theory='dftV3'  %#ok
onlyoriginal=1;  % process db with only original conformations
dbsuffix=''
workname='r9_xyz_dftV3' %#ok %output workname
flexportgeom=1 %#ok choose right output directory
fl_onlyY=0  %export only records with new='Y'
fl_onlyT=0  %export only records with new='T'

savemode.xyz   = 1; 
savemode.zmt   = 0;
savemode.gsxyz = 0; %for this option templatefile is needed
savemode.gszmt = 0; %for this option templatefile is needed
savemode     %#ok

if flexportgeom
    odir=[CD.geomdir filesep workname];
else
    odir=[CD.xyzdir filesep workname];
end

includelist.do=0    %#ok
includelist.type='file';
includelist.name=[odir filesep 'includelist'];
includelist.rectype = 'sdesc' %may have 'sdesc' or 'index' values
includelist  %#ok

excludelist.do=0
%excludelist.type='dir';  % maybe have 'file' or 'dir' values
excludelist.type='file';
excludelist.name=[odir filesep 'excludelist'];
%excludelist.name=[CD.datadir '\g.out\r9\r902.SP\b3lyp'];
excludelist  %#ok

sortmode=0  %#ok %(0/1/2 - without sorting/sort by sdesc/energy)
csconvert=1 %#ok %convert molecules coorginate systems to one
%---------------------------

if ~strcmp(theory,'dft')
  theorystr = ['_' theory];
else                                                                                    
  theorystr = '';
end
workdbname=[CD.dbdir filesep 'r' int2str(moltype) '_g' theorystr];
if onlyoriginal
    templ='_or';
    workdbname = [workdbname templ];
end
if ~isempty(dbsuffix)
    workdbname = [workdbname '_' dbsuffix];
end
workdbname=[workdbname '.mat']  %#ok



%orderanchor=[6 1 2]; %pinds pO4, pC1, pC2

if moltype>=100 && (mod(moltype,100)==11) %tio nucleosides
    orderanchor=[4 70 1]; %pinds pC4, pS4, pC1
else
    orderanchor=[4 6 1]; %pinds pC4, pO4, pC1
end


%orderanchor=[]; %temporarely
construct=0 %#ok
changeisotopes=0;
flPlot=0;
%----------------------------


gtemplname=[workname '_templ.gjf']  %#ok

usedpackage='Gaussian'      %#ok
%usedpackage='xyz'

fullgtemplname=[CD.templatesdir filesep gtemplname] %#ok
if exist(odir,'dir')~=7
    mkdir(odir);
end
tic
load(workdbname,'workdb')

recnum=numel(workdb);
GOenergy=[];
sdesc={};
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

exclset={};
if excludelist.do
if ~isempty(excludelist.name)

  if strcmp(excludelist.type,'dir')
    exclfiles = dir(excludelist.name);
    numfiles = size(exclfiles,1);
  elseif strcmp(excludelist.type,'file')

    [fid,emessage]=fopen(excludelist.name);
    if fid==-1
      error(['No file ' excludelist.name ' with descriptions to exclude found: ' emessage]);
    end

    tline=fgetl(fid);
    i=1;
    while tline~=-1
      exclfiles(i).name=tline;
      tline=fgetl(fid);
      i=i+1;
    end 
    fclose(fid);
    numfiles=i-1;

  end

  if ~numfiles
    disp('No descriptions found to exclude.')
  else
    for i=1:numfiles
      dlm1=strfind(exclfiles(i).name,'_');
      if isempty(dlm1) 
        dlm1=0;
      end   
      dlm2=strfind(exclfiles(i).name,'.');
      if isempty(dlm2), dlm2=size(exclfiles(i).name,2)+1; end;
      exclset(i) = {exclfiles(i).name((dlm1(end)+1):(dlm2-1))}; %complete set of names
    end
  end
end
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

fileind=0;
for i=ind
 
    ms0=workdb(i);
    if fl_onlyY && ~isempty(ms0.new) && ms0.new~='Y' %| 1 %| ms0.new=='P' | ms0.new=='X'
      continue
    end
    if fl_onlyT && ~isempty(ms0.new) && ms0.new~='T' 
      continue
    end
    if strcmpcellar(exclset,ms0.prop.sdesc), continue, end;

    if ~isempty(inclset) && isempty(strcmpcellar(inclset,ms0.prop.sdesc)), continue, end;
    if ~isempty(inclsetind) && ~sum(inclsetind==i), continue, end;



    fileind=fileind+1;
     
    %-----------------------
    if construct
     cla
    %    ms0=zmt2xyz(ms0);
    %    plotmol(ms0,'b')

        ibH41= ms0.ind(find(find(strcmp(pind.labels,'bH41'))==ms0.pind)); 
        ibH42= ms0.ind(find(find(strcmp(pind.labels,'bH42'))==ms0.pind)); 

    if 0 %rb07
        ms0.beta(ibH41)=-ms0.beta(ibH41);
        ms0.beta(ibH42)=-ms0.beta(ibH42);
    %    ms0=zmt2xyz(ms0);
    %    plotmol(ms0,'m')
    end
    if 1 %rb08
        ms0.beta(ibH41)=180*round(ms0.beta(ibH41)/180);
        ms0.beta(ibH42)=180*round(ms0.beta(ibH42)/180);
    %    ms0=zmt2xyz(ms0);
    %    plotmol(ms0,'r')
    end

    %    pause
    end
    if 0
      if isfield(ms0.freq,'freq') && ms0.freq.freq(1)~=0
        [ms0.prop.sdesc ': processing skipped due to freq results already exists']  %#ok
        continue
      end
    end
    %-----------------------
    if changeisotopes % nothing more than change H5 aton to Deuterium
        ibH5= ms0.ind(find(find(strcmp(pind.labels,'bH5'))==ms0.pind)); 
        ms0.masses(ibH5)=2;
    end
    %-----------------------


    if isempty(inclsetind)
        fileindstr=sprintf('%03d',fileind);
    else
        fileindstr=int2str(i);
    end
    ms0.desc=[workname,'#r',int2str(moltype),'_',fileindstr,'_',ms0.prop.sdesc];

    order=[];
    for I=1:numel(orderanchor)
      order(end+1)=find(ms0.pind==orderanchor(I));
    end

%not needed : createzmt now contains createbondchain
%    [xxx,order1]=setdiff(ms0.pind,orderanchor); %output xyz with chosen atoms first then in common way - sorted by pind
%    %ATTENTION: setdiff outputs ROW vector with sorted values
%    order=[order order1];

    %convert all geometries to one coordinate system
    if csconvert
      ms0=createzmt(ms0,order);
      ms0=zmt2xyz(ms0,order);
    end

    if savemode.xyz
      savemol(odir,ms0,0,order);
    end
    if savemode.zmt
      savemol(odir,ms0,3,order);
    end
%    [xxx,order]=sortrows(ms0.pind); %output zmt in common order
    if savemode.gsxyz
      savemolgs(odir,ms0,3,order,fullgtemplname); %Gaussian with XYZ
    end
    if savemode.gszmt
      savemolgs(odir,ms0,5,order,fullgtemplname); %Gaussian with ZMT with varlist
    end

end

toc
