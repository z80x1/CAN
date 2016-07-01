%rot11_3_extrlastgeom: extract last geometry from output gaussian files
%
% Version 1.0    
% Last modified  R O Zhurakivsky 2009-05-21
% Created        R O Zhurakivsky 2009-05-21
time0 = cputime;
tic

format compact
%global flplot

clear ms0 workdb desc new 
atomsind


%--------- data section -------------------------------
workname='rx20' %#ok 
indir='D:\_diplom\CAN2\data\rx19\AIM.2do\' %#ok
odir=[CD.xyzdir filesep workname] %#ok
scanparam='D48';

createdb = 1;
usedpackage='Gaussian'  %#ok
theory='dftV3'          %#ok theory level used for GO
%--------- data section end ---------------------------


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
indir %#ok
diaryfname %#ok

%--------- data section -------------------------------
savemode.xyz   = 0; 
savemode.gsxyz = 0; %for this option templatefile is needed
savemode %#ok
%--------- data section end ---------------------------

gtemplname=[workname '_templ.gjf']  %#ok
fullgtemplname=[CD.templatesdir filesep gtemplname] %#ok


sfiles = dir(strcat(indir,filesep,'*.out'));
sfiles = [sfiles; dir(strcat(indir,filesep,'*.log'))];
numfiles = size(sfiles,1);
if ~numfiles
  error('no files found');
end

for i=1:numfiles

    dlm=strfind(sfiles(i).name,'.');
    fnameshort = sfiles(i).name(1:(dlm-1));
    fnamefull = sfiles(i).name(1:(dlm(end)-1));
    ffname = strcat(indir,filesep,sfiles(i).name); 
    worktitle = lower(fnamefull);
    disp(['Loaded file: ' ffname])

    pat = '_r([a-n\d]+)_';
    A = regexp(fnameshort, pat, 'tokens','once');
    if (numel(A{1})>1) || (A{1}<='9')
        moltype=sscanf(A{1},'%d');
    else
        moltype=10+A{1}-97; %97 - char 'a' code
    end
    pat = '_([A-J][abc]+[ABCSTV])_';
    A = regexp(fnameshort, pat, 'tokens','once');
    conf=A{1};
    
    workdbname=['r' int2str(moltype)];
    workdbname=[workdbname '_' conf]; %#ok
    if strcmp(usedpackage,'Gaussian')
      workdbname=[workdbname '_g']; %#ok
    end
    if ~strcmp(theory,'dft')
      workdbname=[workdbname '_' theory]; %#ok
    end
    workdbname=[CD.dbdir filesep workdbname '.mat']         %#ok

%    try
        fid=fopen(ffname,'r'); % don't use text mode - this is dangerous for Unix files ;)

        frewind(fid);
        workdb=cell(0);
        fl_OptimizedParameters=0;
        fl_OptimizationCompleted=0;
        while 1
            tline = fgetl(fid);
            if ~ischar(tline), fclose(fid); break, end
            if ~isempty(strfind(tline,'Optimization completed'))
              if isempty(ms0)
                  error('Something wrong');
              end

              fl_OptimizationCompleted=1;
            end
            
            if fl_OptimizationCompleted && ~isempty(strfind(tline,'Optimized Parameters'))
                fl_OptimizedParameters=1;
            end
            if ~isempty(strfind(tline,'GradGrad')) %clear all flags
                fl_OptimizedParameters=0;
            end

            if fl_OptimizedParameters && ~isempty(strfind(tline,scanparam))
              pat = '\)\s+(-?\d*\.?\d*)\s+'; %first float number after ')'
              A = regexp(tline, pat, 'tokens','once');
              param=sscanf(A{1},'%f');

              if ~mod(round(param),2) %save only structures with even param value
                  ms0.param=param;

                  if createdb
                      ms0.date=now;
                      ms0 = createbondtable(ms0);
                      [ms0,status] = identmol(ms0,moltype);
                      ms0=calcproperties(ms0,moltype);
                      ms0.comment=param;
                  end
                  if isempty(workdb)
                      workdb=ms0;
                  else
                      workdb(end+1)=ms0; %#ok
                  end
                  disp(param);
              end
            end
            if ~isempty(strfind(tline,'Input orientation'))
              for ii=1:2, tline=fgets(fid); end;  %skip 2 lines to reach "Center     Atomic" line
              [ms0,status]=extrgeom(fid,1,0);
              if status
                  disp([worktitle,': ',lastwarn]);
                  fclose(fid);
                  continue
              end
            end
            if ~isempty(strfind(tline,'SCF Done'))
              [xxx,energy]=strread(tline,'%s%f','delimiter','=');
              ms0.GO.energy=energy;
            end

            
        end

        created=0;
        for num=1:numel(workdb)
            ms0=workdb(num);
            ms0.desc = [fnameshort '_' num2str(ms0.param,'%+04.0f')];

            order=1:ms0.atomnum;
            if savemode.xyz
                savemol(odir,ms0,0,order);
                created=created+1;
            end
            if savemode.gsxyz
                savemolgs(odir,ms0,3,order,fullgtemplname); %Gaussian with XYZ
                created=created+1;
            end

        end
        if created
            disp([created ' gjf files created']);
        end

        if createdb
           dlm=strfind(workdbname,'.');                                          
           workdbnamebkp=[workdbname(1:dlm(end)-1) '~' workdbname(dlm(end):end)];
           if exist(workdbname,'file')                                           
               copyfile(workdbname,workdbnamebkp);                               
           end                                                                   
           save(workdbname,'workdb')                                             
           disp(['Saved db: ' workdbname])
        end
%    catch
%        disp(['catch error: ',lasterr]);
%        fclose(fid);
%    end
end %for

diary off
cputime-time0   %#ok
toc

