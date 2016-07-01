%rot11: analyze output nwchem files
%
% Version 1.0    
% Last modified  R O Zhurakivsky 2007-08-04
% Created        R O Zhurakivsky 2005-09-?

indir=[CD.datadir filesep 'nw.out'];

global pind

clear ms0 dbrib desc
format compact
%grid on
atomsind
maxind=0;

%odir='cyt.out';

load 'dbrib.mat';
%dbrib(1).GO.energy=0;
%dbrib(1).GO.steps=0;
%dbrib(1).SP.energy=0;
%dbrib(1).DM=0;
%dbrib(1).new='';
%dbrib(1).nwfilename='';
%dbrib(1).nwworktitle='';

decssize = numel(dbrib(1).prop.sdesc);
for i=1:numel(dbrib)

%    if isempty(dbrib(i).new),         dbrib(i).new='';, end
%    if isempty(dbrib(i).nwfilename),  dbrib(i).nwfilename='';, end
%    if isempty(dbrib(i).nwworktitle), dbrib(i).nwworktitle='';, end

    desc(i,1:decssize) = dbrib(i).prop.sdesc(1:decssize);
    sdescnew(i) = {dbrib(i).prop.sdescnew};
    nwworktitles(i) = {dbrib(i).nwworktitle};
end

badfilelist=cell(0);

sfiles = dir(strcat(indir,filesep,'*.out'));
numfiles = size(sfiles,1);

maxind=0;
blen=zeros(numfiles,1);
for i=1:numfiles

%    delete(gca);

    dlm=strfind(sfiles(i).name,'.');
    fnameshort = sfiles(i).name(1:(dlm-1));
    fnamefull = sfiles(i).name(1:(dlm(end)-1));

    ffname = strcat(indir,filesep,sfiles(i).name); 

    strcat('Loaded file: ',ffname)

    fid=fopen(ffname,'rt');

    while 1
      tline = fgetl(fid);
      if ~ischar(tline), break, end
      if ~isempty(strfind(tline,'start'))
        break
      end
    end
    [buf,nwworktitle]=strread(tline,'%s%s');

    fseek(fid,0,'eof');
    fseek(fid,-3000,'cof');

    streof = fscanf(fid,'%s',inf);
    if ~isempty(strfind(streof,'CITATION'))

      frewind(fid);

      if strcmpcellar(nwworktitles,nwworktitle)
        fclose(fid);
        continue
      end
      
      while 1
        tline = fgetl(fid);
        if ~ischar(tline), break, end
        if ~isempty(strfind(tline,'Optimization converged'))
           break
        end
      end
      while 1
        tline = fgetl(fid);
        if ~ischar(tline), break, end
        if ~isempty(strfind(tline,'@'))
           break
        end
      end
      buf=sscanf(tline,'@%d%f');
      GO.steps=buf(1);
      GO.energy=buf(2);

      while 1
        tline = fgetl(fid);
        if ~ischar(tline), break, end
        if ~isempty(strfind(tline,'Geometry '))
%           indlastgeom = ftell(fid);
            break
        end
      end
%      fseek(fid,indlastgeom,'bof');
      for i=1:6, fgetl(fid);, end  %skip 6 lines to reach coordinates 
 
      A=fscanf(fid,'%d%s%f%f%f%f',[6,inf]);
%      A = textscan(fid,'%s',10)
      ms0.labels=char(A(2,:))';
      ms0.x=A(4,:)';
      ms0.y=A(5,:)';
      ms0.z=A(6,:)';
      ms0.atomnum = uint16(length(ms0.labels));
      ms0.ind = ((maxind+1):(maxind+ms0.atomnum))';
      ms0.desc = fnameshort;

      while 1
        tline = fgetl(fid);
        if ~ischar(tline), break, end
        if ~isempty(strfind(tline,'Total DFT energy'))
           SP.energy=sscanf(tline,'         Total DFT energy =%f');
           break 
        end
      end
      while 1
        tline = fgetl(fid);
        if ~ischar(tline), break, end
        if ~isempty(strfind(tline,'Dipole moment')) & ...
           ~isempty(strfind(tline,'Debye'))
           DM=sscanf(tline,'   Dipole moment%fDebye(s)');
           break
        end
      end

      ms0 = createbondtable(ms0);
      [ms0,status] = identmol(ms0,,moltype);
      if status
		lastwarn
        fclose(fid);
        continue
      end
      ms0 = createzmt(ms0);
  %    ms0=zmt2xyz(ms0);
      ms0 = calcproperties(ms0,moltype);

      ms0.GO=GO;
      ms0.SP=SP;
      ms0.DM=DM;

      if ms0.prop.hta0<150  %exclude beta-ribose conformations (150 degree is critical value)
        ms0.new='B';
      elseif strcmpar(desc,ms0.prop.sdesc)
        ms0.new='N';
      else
        ms0.new='Y';
        desc(end+1,:)=ms0.prop.sdesc;
        ms0.new
      end

      ms0.nwfilename=fnamefull;
      ms0.nwworktitle=nwworktitle{1};

      ms0.ZPE=[];

      dbrib(end+1)=ms0;

      fclose(fid);
    else %error occured while processing this job
      strcat(nwworktitle,': error')
      badfilelist(end+1)=nwworktitle;

      indlastnonspace=0;
      indlastspace=0;
      indlastspace2=0;
      fseek(fid,-3000,'eof');
      while 1
        tline = fgetl(fid);
        if ~ischar(tline), break, end
        if ~numel(tline) | tline(1)==' '
          indlastspace2 = indlastspace;  
          indlastspace = ftell(fid);
        else
          indlastnonspace = ftell(fid);
        end

      end
      if indlastspace>indlastnonspace
        indlastspace=indlastspace2;
      end
      fseek(fid,indlastspace,'bof');
      tline = fgetl(fid);
      tline = fgetl(fid);

      A=find(tline==':');
      B=circshift(A,[0 -1]);
      C=B-A;
      ind=find(C>5);
      if isempty(ind)
          erdirname='unknown';
      else
          erdirname=tline(A(ind(1))+1:A(ind(1)+1)-1);
      end
      erdirname=strcat(indir,filesep,erdirname);
      if ~isdir(erdirname)
        mkdir(erdirname);
      end

      fclose(fid);
      movefile(ffname,erdirname,'f');
%      movefile(ffname,strcat(erdirname,filesep,sfiles(i).name));

    end
 
%    'Press any key.'
%    pause
end

save 'dbrib.mat' dbrib
