%rot11: analyze output gaussian files
%
% Version 1.0    
% Last modified  R O Zhurakivsky 2008-12-21
% Created        R O Zhurakivsky 2005-09-?

%070501: L528: changed errorous i index to n
%090915: added moltype == 910 support

clear
format compact
global pind GLaspec
global gl_fl_oldchistyle

%clear ms0 workdb desc new freq ZPE
pindsdef
atomsind
dbchanged=0;
dbchangedunique=0;

%--------- data section -------------------------------
%CD.dbdir='D:\zhr\_diplom\work\Alla\090114\r312';
%CD.datadir=CD.dbdir;
%--------- data section end ---------------------------

%indir=[CD.datadir filesep 'r910/r91001'];
indir='D:\_diplom\temp\r91006_restart_bad_sugar.gator';

diaryfname0=[indir filesep 'logfile'];
diaryfname=diaryfname0;
for i=2:inf
  if ~(exist(diaryfname,'file')==2)
    break
  end
  diaryfname = [diaryfname0 int2str(i)];
end
diary(diaryfname)


disp(['Time: ' datestr(now)]);
indir                   %#ok
diaryfname              %#ok

%--------- data section -------------------------------
moltype=910             %#ok
usedpackage='Gaussian'  %#ok
%theory='cryst'         %#ok 
theory='dftV3'          %#ok theory level used for GO
%dftV2 - b3lyp/6-31G(d,p)
%dftV3 - b3lyp/6-31G(d,p) DB with 4 syn/anti states
%dftV4 - b3lyp/6-31G++(d,p)
%dftV5 - b3lyp/6-311G++(d,p)
%mp2V2 - mp2/6-31G(d,p)
%090506 r11_g_dftV3.mat renamed to r11_g_dftV2.mat
dbsuffix=''

flexacttitle=1          %use for work identification exact match of worknames or inclusion only

updatedb      = 0      %#ok %true if updating existing structures
savedb        = 1       %#ok %true for save database after processing
processerrors = 0       %#ok %true for processing datafiles abnormally terminated
task.opt      = 1;      %true if files to process consist GO data
task.freq     = 0;      %true if files to process consist FC data
task.energy   = 0;      %true if files to process consist SP data
%task.energyfieldname='MP2_6311Gdp';
task.energyfieldname='MP2_6311__Gdp';
%task.energyfieldname='MP2_6311__G2dfpd'; %name of field to write energy data in
%task.energyfieldname='B3LYP_6311Gdp';
flfreqchk       = 0       %true if using freqchk utility
task.energyMP2  = 1;      %true if SP data is MP2 energy data
task.updateTS   = 0;      %true if files to process consist data on trasition state structures
task.updateFixed= 0;      %true if files to process consist data on trasition state structures ???
task.updateDNA  = 0;      %update records with new='D'
task                    %#ok

gl_fl_oldchistyle= 0; %if is set is 'A' symbol for all anti- conformers, and 'S' for all syn-
%--------- data section end ---------------------------


if ~task.freq 
    flfreqchk=0;
end

if moltype<30
      workdbname=['r' int2str(moltype)];
elseif moltype>=30 && moltype<100
      workdbname=['r' sprintf('%03.0f',moltype)];
elseif moltype>=100
      workdbname=['r' int2str(moltype)];
else
  warning('rot11:unktype','unknown molecule type')
  break
end

if strcmp(usedpackage,'Gaussian')
  workdbname=[workdbname '_g'];
end
if ~strcmp(theory,'dft')
  workdbname=[workdbname '_' theory];
end

newworkdbname=[workdbname '_or'];   %for DB with original conformations
if ~isempty(dbsuffix)
    newworkdbname = [newworkdbname '_' dbsuffix];
    workdbname = [workdbname '_' dbsuffix];
end

newworkdbname=[CD.dbdir filesep newworkdbname '.mat']   %#ok
workdbname=[CD.dbdir filesep workdbname '.mat']         %#ok


if exist(workdbname,'file')==2
  load(workdbname,'workdb');
else
  warning('rot11:dbnexist','Database doesn''t exists!')  
  workdb={};  
end

sdescnew={}; worktitles={}; new=[];
desc=char();
if ~isempty(workdb)
  decssize = numel(workdb(1).prop.sdesc);
  recnum=numel(workdb);
  for i=1:recnum
    desc(i,1:decssize) = workdb(i).prop.sdesc(1:decssize); %conformation identification
    if isfield(workdb(i).prop,'sdescnew')
        sdescnew(i) = {workdb(i).prop.sdescnew};
    end
    worktitles(i) = {lower(workdb(i).worktitle)};
    new(i) = workdb(i).new;
  end
else
%  worktitles={};
  new=zeros(1,0);
  sdescnew={};
  worktitles={};
%  workdb(1).worktitle='';
end

badfilelist=cell(0);

sfiles = dir(strcat(indir,filesep,'*.out'));
sfiles = [sfiles; dir(strcat(indir,filesep,'*.log'))];
numfiles = size(sfiles,1);
if ~numfiles
  error('No files found');
end

blen=zeros(numfiles,1);
for i=1:numfiles

%    delete(gca);
    ms0=struct();

    dlm=strfind(sfiles(i).name,'.');

    fnameshort = sfiles(i).name(1:(dlm-1));
    dlm2=strfind(fnameshort,'_');
    if ~isempty(dlm2)
        fnameshort=fnameshort(dlm2(end)+1:end);
    end

    fnamefull = sfiles(i).name(1:(dlm(end)-1));

    ffname = strcat(indir,filesep,sfiles(i).name); 

    disp(['Loaded file: ' ffname])

    worktitle = lower(fnamefull);
    worktitleshort=lower(fnameshort);
    fid=fopen(ffname,'r'); % don't use text mode - this is dangerous for Unix files ;)

    fseek(fid,-500,'eof');
    streof = fscanf(fid,'%s',inf);

    erroroccured=isempty(strfind(streof,'Normaltermination'));
    if processerrors || ~erroroccured || flfreqchk

      if flfreqchk
        comment = 'updated with flfreqchk data';    %#ok
      elseif erroroccured 
        comment = 'Error during calculation'    %#ok
      else
        comment='';
      end

      if flexacttitle
          indwork = strcmpcellar(worktitles,worktitle,flexacttitle); %index of work with same name
      else % try to find similar work
          indwork = strcmpcellar(worktitles,worktitleshort,flexacttitle); %index of work with same name
      end
      
      if ~(isempty(indwork) || updatedb) %work with the same name exists, don't add this one to end of base 
          disp(['warning: ' worktitle ' exists. continue'])
          fclose(fid);
          continue
      else
          dbchanged=1; %if we are here something must be added to db
      end

      if updatedb % if updating base

         dbchangedunique=1; %if db is updated there is large probability of uniqueconfdb to reshesh 
         if 0
           if isempty(indwork) % no work with the same name, maybe one with same desc field? test this
             dlm=strfind(fnameshort,'_');
  %           dlmminus=strfind(fnameshort,'-');
             if isempty(dlm)
               dlm=0;
             end
  %           dlm=strfind(sfiles(i).name,'_');
  %          if dlmminus(end) > dlm(end)
  %            iend = dlmminus(end)-1;
  %          else
               iend = numel(fnameshort);
  %          end
             thisdesc = fnameshort((dlm(end)+1):iend);
             buf=repmat(thisdesc,size(desc,1),1);
             indwork = find(sum(buf==desc,2)==numel(thisdesc) & new'=='Y');
           end                             
         end

         frewind(fid);
         if ~flfreqchk
           [ms0,status]=extrgeom(fid,0,1);
         else %processing freqchk datafile
           for ijk=1:2, tline=fgets(fid); end  %skip 2 lines to reach "Center     Atomic" line
           [ms0,status]=extrgeom(fid,moltype,1);
         end
         if status
            disp([worktitle,': ',lastwarn])
            fclose(fid);
            continue
         end
         ms0 = createbondtable(ms0);
         [ms0,status] = identmol(ms0,moltype);
         if status
            disp([worktitle,': ',lastwarn]);
            fclose(fid);
            continue
         end

         % ms0 = createzmt(ms0);
         % ms0=zmt2xyz(ms0);

         if moltype == 910
             ms0 = calcproperties910(ms0,moltype);
         else
             ms0 = calcproperties(ms0,moltype);
         end

         ms0.GO.steps=0;
         ms0.SP=0;
         %ms0.gaussian=gaussian;

         thisdesc=ms0.prop.sdesc;
         if size(desc,1)>0
              buf=repmat(thisdesc,size(desc,1),1);
              if task.updateTS|task.updateFixed
                  if task.updateTS & isempty(indwork)
                      indwork = find(sum(buf==desc,2)==numel(thisdesc) & (new'=='T') );
                  end
                  if task.updateFixed & isempty(indwork)
                      indwork = find(sum(buf==desc,2)==numel(thisdesc) & (new'=='F')|(new'=='C') );
                  end
              elseif task.updateDNA
                  if isempty(indwork)
                      indwork = find(sum(buf==desc,2)==numel(thisdesc) & (new'=='D') );
                  end
              else
                  indwork = find(sum(buf==desc,2)==numel(thisdesc) & new'=='Y');
              end
         end
         if isempty(indwork)% no similar conformations in the base
            disp(['error01: ' thisdesc ' not found for updating'])
            fclose(fid);
            continue
         end

         dx=0.0001;
%        dr=0.0008; da=0.01; dt=0.01;
         dr=0.0008; da=0.1; dt=1;

         ms0=createzmt(ms0,find(ms0.pind==strcmpcellar(pind.labels,'pC1')));

         indworkgood=[];
         maxdall=[];
         for indwork_=reshape(indwork,1,numel(indwork))

%            indwork=0;
             flmolidentical=1;
             mstmp=createzmt(workdb(indwork_),find(workdb(indwork_).pind==strcmpcellar(pind.labels,'pC1')));
        
               [xxx,mstmpI]=sort(mstmp.pind);  %#ok
               [xxx,ms0I]=sort(ms0.pind);

               if ~(mstmp.iR(mstmpI)==ms0.iR(ms0I) &...
                 mstmp.ialfa(mstmpI)==ms0.ialfa(ms0I) &...
                 mstmp.ibeta(mstmpI)==ms0.ibeta(ms0I))%#ok
                 flmolidentical=0;
               end
               drcur = abs(mstmp.R(mstmpI)-ms0.R(ms0I));
               if sum(drcur>dr) %sum is for making OR
                 flmolidentical=0;
               end
               dacur = abs(mstmp.alfa(mstmpI)-ms0.alfa(ms0I));
               if sum(dacur>da)
                 flmolidentical=0;
               end
               dtcur = abs(mstmp.beta(mstmpI)-ms0.beta(ms0I));
               dtcur=abs(dtcur-floor((dtcur+dt)/360)*360); %avoiding situation when one angle is -179.9999 and second is 180.0000
               if sum(dtcur>dt)
                  flmolidentical=0;
               end

             if flmolidentical==1
                 indworkgood(end+1) = indwork_;
                 maxdall(end+1) = max(drcur)*max(dacur)*max(dtcur);
                 break;
             end
        
         end
         
         if numel(indworkgood)
            [Y,I]=min(maxdall);
            indwork=indworkgood(I);
         else
            disp([worktitle ' has conformation not identical to ' thisdesc ' (#' int2str(indwork) ') ['...
                  num2str(max(drcur),2) ',' num2str(max(dacur),2) ',' num2str(max(dtcur),2) ']'])
            indwork=0;
         end
        if ~indwork % no similar conformations in the base
          disp(['error02: ' thisdesc ' not found for updating'])
          fclose(fid);
          continue
        end

        ms0=workdb(indwork);
        comment=['Updated with ' worktitle ' data at ' datestr(now)];
        disp([worktitle ' (#' int2str(indwork) ',' ms0.prop.sdesc ') : updating... [' ...
              num2str(max(drcur),2) ',' num2str(max(dacur),2) ',' num2str(max(dtcur),2) ']'])

      end % if updatedb

      frewind(fid);
      posStandGeom=0;
     
      if task.opt %analyse geometry optimization results

          if updatedb
            gaussian=ms0.gaussian;
          end

          while 1
            tline = fgetl(fid);
            if ~ischar(tline), break, end

            if ~isempty(strfind(tline,'SCF Done'))
              [xxx,B3LYP_631Gdp]=strread(tline,'%s%f','delimiter','=');
            end
%            if ~isempty(strfind(tline,'Standard orientation')) %zhr090522
            if ~isempty(strfind(tline,'Input orientation'))
              for ijk=1:2, tline=fgets(fid); end  %skip 2 lines to reach "Center     Atomic" line
              posStandGeom=ftell(fid);
            end
%            if ~isempty(strfind(tline,'Optimization completed'))
            if ~isempty(strfind(tline,'Optimization '))
                break
            end
          end
          tline = fgetl(fid);
          if ~isempty(strfind(tline,'Stationary point found'))
            stpoint='Y';
          else
            stpoint='N';
          end

          fseek(fid,posStandGeom,'bof');
          [ms0,status]=extrgeom(fid,moltype,1);
          if status
            disp([worktitle,': ',lastwarn]);
            fclose(fid);
            continue
          end

          ms0 = createbondtable(ms0);
          [ms0,status] = identmol(ms0,moltype);
          if status
            disp([worktitle,': ',lastwarn]);
            fclose(fid);
            continue
          end

          ms0 = createzmt(ms0);
          % ms0=zmt2xyz(ms0);
          if moltype == 910
              ms0 = calcproperties910(ms0,moltype);
          else
              ms0 = calcproperties(ms0,moltype);
          end

          ms0.GO.steps=0;
          ms0.SP=0;
          %ms0.gaussian=gaussian;

          ms0.GO.energy=B3LYP_631Gdp;
          ms0.gaussian.B3LYP_631Gdp=B3LYP_631Gdp;
          ms0.gaussian.stpoint=stpoint;

          if ~updatedb
              [ms0,desc]=checknewconf(ms0,workdb,desc);
              if ms0.new=='Y'
                  dbchangedunique=1;
              end
              ms0.desc = fnameshort;
              ms0.filename=fnamefull;
              ms0.worktitle=worktitle;
          end

 
      end

      if task.energy
%keyboard
          if isempty(fieldnames(ms0))
            disp('No molecule structure loaded, skipping')
            fclose(fid);
            continue
          end

          if ~task.energyMP2

            while 1
              tline = fgetl(fid);
              if ~ischar(tline)
                disp([worktitle,': B3LYP energy not updated.'])
                break
              end
              if ~isempty(strfind(tline,'SCF Done'))
                  [xxx,energy]=strread(tline,'%s%f','delimiter','=');
                  ms0.gaussian.(task.energyfieldname)=energy;
                  break;
              end
            end  

          else %MP2

            while 1
              tline = fgetl(fid);
              if ~ischar(tline)
                disp([worktitle,': MP2 energy not updated.'])
                break
              end
              if ~isempty(strfind(tline,'EUMP2'))
                  [xxx,xxx2,xxx3,energy]=strread(tline,'%s%f%s%f','delimiter','=');
                  ms0.gaussian.(task.energyfieldname)=energy;
                  break;
              end
            end

          end

          while 1
            tline = fgetl(fid);
            if ~ischar(tline)
              disp([worktitle,': Charges information not found.'])
              break
            end
            if ~isempty(strfind(tline,'Mulliken atomic charges:'))
                tline = fgetl(fid);  %skip 1 lines to reach charges 

                A=fscanf(fid,'%d %s %f',[3,inf]);
                if isempty(A)
                  warning('rot11:outputerror','fatal error in output file - couldn''n load charges. skipping');
                  break
                end
%                ms0.gaussian=setfield(ms0.gaussian,'mcharge',A(3,:)); %Mulliken atomic charges
                ms0.gaussian.('mcharge')=A(3,:); %Mulliken atomic charges

                break;
            end
          end

      end

      if task.freq %analyse frequency results

          if isempty(fieldnames(ms0))
            disp('No molecule structure loaded, skipping')
            fclose(fid);
            continue
          end

          frewind(fid);
          freq=[];
%keyboard

          while 1
            tline = fgetl(fid);
            if ~ischar(tline), break, end
            flgeomfound=0;
            if flfreqchk 
              if ~isempty(strfind(tline,'Center     Atomic'))
                flgeomfound=2;
              end
            else
%              if ~isempty(strfind(tline,'Standard orientation')) %zhr090522
              if ~isempty(strfind(tline,'Input orientation')) 
                flgeomfound=1;
              end
            end

            if flgeomfound
                if flgeomfound==1 %standart work mode
                    for ijk=1:2, fgetl(fid); end  %skip 2 lines to reach "Center     Atomic" line
                end
                tline=fgets(fid); %now we are at line "Number     Number                        X           Y           Z"
                numfields  = numel(strread(tline,'%s'));
                tline=fgets(fid); %skip line to reach coordinates 

                if numfields==5
                    A=fscanf(fid,'%d%d%f%f%f',[5,inf]);
                elseif numfields==6
                    A=fscanf(fid,'%d%d%d%f%f%f',[6,inf]);
                else
                    error('incorrect number of fields in geometry structure detected');
                end
                if isempty(A)
                  warning('rot11:outpuiterror2','fatal error in output file - couldn''t load freq geometry. skipping');
                  continue
                end
                [a,b]=meshgrid(A(2,:),GLaspec.weight);
                c=a==b;
                for ijk=1:numel(GLaspec.weight)
                  freq.labels(c(ijk,:))=GLaspec.type(ijk);
                end
                if ~isempty(strcmpcellar(freq.labels,''))
                  warning('rot11:emptylabels','Empty labels found!');
                end
                freq.x=A(numfields-2,:)';
                freq.y=A(numfields-1,:)';
                freq.z=A(numfields,:)';
                break
            end
          end  

          while 1
            tline = fgetl(fid);
            if ~ischar(tline), break, end
            if ~isempty(strfind(tline,'Harmonic frequencies')) 
               
               n=1;
               while n < ms0.atomnum*3-6
                 while 1
                   if ~ischar(tline), break, end
                    if strncmp(tline,'               ',15) 
                      break
                    end
                    tline = fgetl(fid);
                 end
                 tline = fgetl(fid); %#ok
                 tline = fgetl(fid); 
                 [buf,freq.freq(n),freq.freq(n+1),freq.freq(n+2)]=strread(tline,'%16c%f%f%f');
                 tline = fgetl(fid); %#ok
                 tline = fgetl(fid); 
                 [buf,freq.forceconst(n),freq.forceconst(n+1),freq.forceconst(n+2)]=strread(tline,'%16c%f%f%f');%#ok
                 tline = fgetl(fid);
                 [buf,freq.intencity(n),freq.intencity(n+1),freq.intencity(n+2)]=strread(tline,'%16c%f%f%f');

                 while 1
                   tline = fgetl(fid);
                   if ~ischar(tline), break, end
                   if ~isempty(strfind(tline,'Atom AN'))
                      tline = fgetl(fid);
                      break
                   end
                 end
                 for j=1:ms0.atomnum
                   [xxx,xxx2,freq.dx(n,j),freq.dy(n,j),freq.dz(n,j),freq.dx(n+1,j),freq.dy(n+1,j),freq.dz(n+1,j),freq.dx(n+2,j),freq.dy(n+2,j),freq.dz(n+2,j)]=strread(tline,'%d%d%f%f%f%f%f%f%f%f%f');
                   tline = fgetl(fid);
                 end
                 n=n+3;

              end
               break
            end
          end

          [ZPE,TEC,EEC,GEC,IA,IB,IC]=deal(inf);
          while 1
            tline = fgetl(fid);
            if ~ischar(tline), break, end

            if ~isempty(strfind(tline,'Temperature')) 
               [xxx,T]=strread(tline,'%s%f',1,'delimiter',' '); %Kelvin
            end

            if ~isempty(strfind(tline,'Principal axes and moments of inertia')) 
            tline = fgetl(fid);%#ok
            tline = fgetl(fid);
            [buf]=strread(tline,'%s');
                if any(buf{3}=='*')
                    comment=strcat(comment,', Moments of inertia are over limits!');
                    disp('Warning: moments of inertia are over limits!');
                    [IA,IB,IC]=deal(inf);
                else
                    [IA,IB,IC]=strread(buf{3},'%10f%10f%10f');
                end
            end
            if ~isempty(strfind(tline,'Zero-point correction')) 
               [xxx,ZPE]=strread(tline,'%s%f','delimiter','=');
               ZPE=ZPE*CC.encoef;  %kcal/mol
            end
            if ~isempty(strfind(tline,'Thermal correction to Energy')) 
               [xxx,TEC]=strread(tline,'%s%f','delimiter','=');
               TEC=TEC*CC.encoef;  %kcal/mol
            end
            if ~isempty(strfind(tline,'Thermal correction to Enthalpy')) 
               [xxx,EEC]=strread(tline,'%s%f','delimiter','='); %Enthalpy energy correction
               EEC=EEC*CC.encoef;  %kcal/mol
            end
            if ~isempty(strfind(tline,'Thermal correction to Gibbs Free Energy')) 
               [xxx,GEC]=strread(tline,'%s%f','delimiter','='); %Gibbs energy correction
               GEC=GEC*CC.encoef;  %kcal/mol 
               break
            end

          end

          ms0.ZPE=ZPE;
          ms0.gaussian.ZPE=ZPE; %nonscaled values

          if (isfield(ms0.gaussian,'TEC') || isfield(ms0.gaussian,'EEC') || isfield(ms0.gaussian,'GEC')) && ~isfield(ms0.gaussian,'T')
            ms0.gaussian.T=298.15;
          end
          
          if ~isfield(ms0.gaussian,'T')
            ms0.gaussian.T=T;
            ms0.gaussian.TEC=TEC;
            ms0.gaussian.EEC=EEC;
            ms0.gaussian.GEC=GEC;
          else
            tind = find(ms0.gaussian.T==T);
            if tind
                ms0.gaussian.T(tind)=T;
                ms0.gaussian.TEC(tind)=TEC;
                ms0.gaussian.EEC(tind)=EEC;
                ms0.gaussian.GEC(tind)=GEC;
            else
                ms0.gaussian.T(end+1)=T;
                ms0.gaussian.TEC(end+1)=TEC;
                ms0.gaussian.EEC(end+1)=EEC;
                ms0.gaussian.GEC(end+1)=GEC;
            end
          end

          ms0.gaussian.imoments=[IA IB IC]; %moments of inertia
          ms0.freq=freq;

      elseif updatedb==0
          ms0.ZPE=inf;
          ms0.gaussian.ZPE=inf;

          ms0.gaussian.TEC=inf;
          ms0.gaussian.EEC=inf;
          ms0.gaussian.GEC=inf;

          ms0.gaussian.imoments=[0 0 0]; %moments of inertia
          ms0.freq.freq=0;
      end %freq

      if ~updatedb
        ms0.comment={comment};
        ms0.date=now;
      elseif ~isempty(comment)
        if isempty(ms0.comment(end))
          ms0.comment{end}=comment;
        else
          ms0.comment{end+1}=comment;
        end
      end


      if isempty(workdb)
        workdb=ms0;
        new(1)=ms0.new;
        sdescnew(1) = {ms0.prop.sdescnew};               
        worktitles(1) = {lower(ms0.worktitle)};          
      else
        ms0=orderfields(ms0,workdb(1));

        if updatedb
           workdb(indwork)=ms0;
           disp(['workdb[' int2str(indwork) ']'])

        else
           indwork=numel(workdb)+1;
           workdb(indwork)=ms0;
           disp(['workdb[' int2str(indwork) ']'])

           new(end+1)=ms0.new;
           %desc(end+1,1:decssize) = ms0.prop.sdesc(1:decssize);  %is added by checknewconf
           if isfield(ms0.prop,'sdescnew')
               sdescnew(end+1) = {ms0.prop.sdescnew};               
           else
               sdescnew(end+1) = {''};
           end
           worktitles(end+1) = {lower(ms0.worktitle)};          
        end  
      end
      fclose(fid);

    else %error occured while processing this job
      strcat(worktitle,': output file error')
      badfilelist(end+1)={worktitle};

      fclose(fid);
      erdirname=strcat(indir,filesep,'error');
      if ~isdir(erdirname)
        mkdir(erdirname);
      end
      movefile(ffname,erdirname,'f');
    end
 
%break
%    'Press any key.'
%    pause
end

fl_tosave=0;
fl_ressaved=0;

if savedb && dbchanged 
  if processerrors
    disp('You''ve specified to save non successfully completed file''s data in database');
    answer=input('Are you sure you want to proceed? 1/[0]');
    if answer==1
      fl_tosave=1;
    end
  else
    fl_tosave=1;
  end
end

if fl_tosave
    dlm=strfind(workdbname,'.');
    workdbnamebkp=[workdbname(1:dlm(end)-1) '~' workdbname(dlm(end):end)];
    if exist(workdbname,'file')
        copyfile(workdbname,workdbnamebkp);
    end
    save(workdbname,'workdb')

    
    disp([' Total number of conformations in DB: ' int2str(numel(workdb))])
    disp([' Total number of unique conformations: ' int2str(size(unique(desc,'rows'),1))])
    disp([' Total number of Y conformations: ' int2str(sum(new=='Y'))])

    if dbchangedunique
      disp('Saving db with unique conformations...')
      workdb=workdb(new=='Y');
      dlm=strfind(newworkdbname,'.');
      newworkdbnamebkp=[newworkdbname(1:dlm(end)-1) '~' newworkdbname(dlm(end):end)];
      if exist(newworkdbname,'file')
        copyfile(newworkdbname,newworkdbnamebkp);
      end
      save(newworkdbname,'workdb')
    end

    fl_ressaved=1;
end

if ~fl_ressaved
    disp('Results are not saved!')
else
    clear workdb  
end


clear A a b c a1 a3 a5 a7 xxx 
diary off

 
