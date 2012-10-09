%rot11_4_dbextrfromscan: extract geometries from scan E(chi) output files
%
% Version 1.0    
% Last modified  R O Zhurakivsky 2010-01-11
% Created        R O Zhurakivsky 2009-05-21
time0 = cputime;
tic

format compact
%global flplot

clear ms0 workdb desc new 
atomsind


%--------- data section -------------------------------
workname='rx21_mp2V2_WFN' %#ok 
indir='D:\_diplom\CAN2\data\rx21_mp2V2\r9' %#ok
odir=[CD.xyzdir filesep workname] %#ok
scanparam='D34'; % !!may be found from 'Initial Parameters' section


createdb = 1;
fl_reusedb = 1;
usedpackage='Gaussian'  %#ok
theory='mp2V2'          %#ok theory level used for GO
dbsuffix='chilimited'

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
savemode.gsxyz = 1; %for this option templatefile is needed
savemode %#ok

limits4process = [
%{'rx21_r9_AabcA_tchi_scan_left.gjf.01-j1.out'}  {57}    {360-163} ;
%{'rx21_r9_AabcA_tchi_scan_right.gjf.01-j2.out'} {-162}  {-107};
%{'rx21_r9_EabcA_tchi_scan_left.gjf.01-j3.out'}  {111}   {360-146};
%{'rx21_r9_EabcA_tchi_scan_right.gjf.01-j4.out'} {-145}  {-37};
%{'rx21_rd_AabcA_tchi_scan_left.gjf.01-j5.out'}  {42}    {360-159};
%{'rx21_rd_AabcA_tchi_scan_right.gjf.01-j6.out'} {-158}  {-104};
%{'rx21_rd_EabcA_tchi_scan_left.gjf.01-j7.out'}  {99}    {360-129};
%{'rx21_rd_EabcA_tchi_scan_right.gjf.01-j8.out'} {-128}  {-24};
%{'rx21_rf_AabcA_tchi_scan_left.gjf.01-j9.out'}  {44}    {360-153};
%{'rx21_rf_AabcA_tchi_scan_right.gjf.01-j10.out'} {-152} {-122};
%{'rx21_rf_EabcA_tchi_scan_left.gjf.01-j11.out'}  {106}  {360-125};
%{'rx21_rf_EabcA_tchi_scan_right.gjf.01-j12.out'} {-124} {-33};
%{'rx21_rg_AabcA_tchi_scan_left.gjf.01-j13.out'}  {150}  {360-149};
%{'rx21_rg_AabcA_tchi_scan_right.gjf.01-j14.out'} {-148} {-132}; %zhr100114 max value corrected
%{'rx21_rg_EabcA_tchi_scan_left.gjf.01-j15.out'}  {115}  {360-123};
%{'rx21_rg_EabcA_tchi_scan_right.gjf.01-j16.out'} {-122} {-71};

%commented at 2010-0613
%{'rx21_r9_AabcA_tchi_scan_left.gjf.01-j1.out'}  {179}    {360-163} ;
%{'rx21_r9_AabcA_tchi_scan_right.gjf.01-j2.out'} {-162}  {-139};
%{'rx21_r9_EabcA_tchi_scan_left.gjf.01-j3.out'}  {1000}   {0};
%{'rx21_r9_EabcA_tchi_scan_right.gjf.01-j4.out'} {234}  {272};
%{'rx21_rd_AabcA_tchi_scan_left.gjf.01-j5.out'}  {179}    {360-159};
%{'rx21_rd_AabcA_tchi_scan_right.gjf.01-j6.out'} {-158}  {-139};
%{'rx21_rd_EabcA_tchi_scan_left.gjf.01-j7.out'}  {1000}    {0};
%{'rx21_rd_EabcA_tchi_scan_right.gjf.01-j8.out'} {234}  {272};
%{'rx21_rf_AabcA_tchi_scan_left.gjf.01-j9.out'}  {179}    {360-153};
%{'rx21_rf_AabcA_tchi_scan_right.gjf.01-j10.out'} {-152} {-139};
%{'rx21_rf_EabcA_tchi_scan_left.gjf.01-j11.out'}  {234}  {360-125};
%{'rx21_rf_EabcA_tchi_scan_right.gjf.01-j12.out'} {360-124} {272};
%{'rx21_rg_AabcA_tchi_scan_left.gjf.01-j13.out'}  {179}  {360-149};
%{'rx21_rg_AabcA_tchi_scan_right.gjf.01-j14.out'} {-148} {-139};
%{'rx21_rg_EabcA_tchi_scan_left.gjf.01-j15.out'}  {234}  {360-123};
%{'rx21_rg_EabcA_tchi_scan_right.gjf.01-j16.out'} {-122} {272};
];
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

    fname = sfiles(i).name;
    dlm=strfind(fname,'.');
    fnameshort = fname(1:(dlm-1));
    fnamefull = fname(1:(dlm(end)-1));
    ffname = strcat(indir,filesep,fname); 
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
    if ~isempty(dbsuffix)
        workdbname = [workdbname '_' dbsuffix];
    end

    workdbname=[CD.dbdir filesep workdbname '.mat']         %#ok

%    try
        fl_param_limited = 0;
        for ii=1:size(limits4process,1)
            if  strcmp(limits4process{ii,1},fname)
                param_min = limits4process{ii,2};
                param_max = limits4process{ii,3};
                fl_param_limited=1;
                break;
            end
        end
        fid=fopen(ffname,'r'); % don't use text mode - this is dangerous for Unix files ;)

        frewind(fid);
        
        if createdb && fl_reusedb
            if exist(workdbname,'file')==2
              load(workdbname,'workdb');
            else
              warning('rot11:dbnexist','Database doesn''t exists!')  
              workdb=cell(0);  
            end
        else
            workdb=cell(0);
        end
        
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
              param=sscanf(A{1},'%f'); %scanned parameter's value

              if fl_param_limited
                  if round(mod(param,360))<round(mod(param_min,360))
                      disp(['param value ' int2str(param) ' skipped as less then min']);
                      continue
                  end
                  if round(mod(param,360))>round(mod(param_max,360))
                      disp(['param value ' int2str(param) ' skipped as more then max']);
                      continue
                  end
              end
              
              if 0 && mod(round(param),2) %save only structures with even param value
                  continue;
              end
              
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
            if exist(odir,'dir')~=7
                mkdir(odir);
            end
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
            disp([int2str(created) ' gjf files created']);
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

