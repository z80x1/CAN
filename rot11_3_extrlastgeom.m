%rot11_3_extrlastgeom: extract last geometry from output gaussian files
%
% Version 1.0    
% Last modified  R O Zhurakivsky 2008-06-37
% Created        R O Zhurakivsky 2005-09-?

format compact
global flplot
%global pind GLaspec

clear ms0 workdb desc new 
atomsind
indir='';


molecule0 = 'r312' %#ok
workname0 = 'r31205' %#ok
add_dir = 'error' %#ok

%indir='D:\zhr\_diplom\work\Alla\090110\' %#ok

if ~length(indir)
  indir=[CD.datadir filesep molecule0 filesep workname0 ] %#ok
end

if length(add_dir), indir=[indir filesep add_dir];, end

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
moltype=413 %#ok
usedpackage='Gaussian' %#ok
theory='dftV3' %#ok


%workname=workname0 %#ok
workname='r11205v2'

savemode.xyz   = 0;
savemode.zmt   = 0;
savemode.gsxyz = 1; %for this option templatefile is needed
savemode.gszmt = 0; %for this option templatefile is needed
savemode     %#ok

fl_identmol=0  %#ok %identify or not identify - that is the question
%orderanchor=[6 1 2]; %pinds pO4, pC1, pC2
orderanchor=[4 6 1]; %pinds pC4, pO4, pC1
fl_createmovie=0 %#ok
fl_extractall=0 %#ok %extract geometry at all steps
flplot=1    %#ok

gtemplname=[workname '_templ.gjf']  %#ok
fullgtemplname=[CD.templatesdir filesep gtemplname] %#ok

%-----------------------------------------
workdbname=['r' int2str(moltype)];
if strcmp(usedpackage,'Gaussian')
  workdbname=[workdbname '_g'];
end
if ~strcmp(theory,'dft')
  workdbname=[workdbname '_' theory];
end

workdbname=[CD.dbdir filesep workdbname '.mat'] %#ok

odir=indir;
%odir=[CD.xyzdir filesep workname];
if exist(odir,'dir')~=7
   mkdir(odir);
end

if exist(workdbname,'file')==2
  load(workdbname,'workdb');
else
  workdb={};  
end

if ~isempty(workdb)
  desc=char();
  decssize = numel(workdb(1).prop.sdesc);
  for i=1:numel(workdb)
    desc(i,1:decssize) = workdb(i).prop.sdesc(1:decssize); %conformation identification
    new(i) = workdb(i).new;
  end
else
  desc={};
  new={};
end

sfiles = dir(strcat(indir,filesep,'*.out'));
sfiles = [sfiles; dir(strcat(indir,filesep,'*.log'))];
numfiles = size(sfiles,1);
if ~numfiles
  error('no files found');
end

for i=1:numfiles

%    ms0=struct();

    dlm=strfind(sfiles(i).name,'.');
    fnameshort = sfiles(i).name(1:(dlm-1));
    fnamefull = sfiles(i).name(1:(dlm(end)-1));

    ffname = strcat(indir,filesep,sfiles(i).name); 

    disp(['Loaded file: ' ffname])

    worktitle = lower(fnamefull);


%    try
        fid=fopen(ffname,'r'); % don't use text mode - this is dangerous for Unix files ;)

        frewind(fid);
        posInputGeom=[];
        posStandGeom=[];
        while 1
            tline = fgetl(fid);
            if ~ischar(tline), break, end

            if ~isempty(strfind(tline,'Standard orientation')) 
              for i=1:2, tline=fgets(fid); end  %skip 2 lines to reach "Center     Atomic" line
              posStandGeom(end+1)=ftell(fid);
            elseif ~isempty(strfind(tline,'Input orientation')) 
              for i=1:2, tline=fgets(fid); end  %skip 2 lines to reach "Center     Atomic" line
              posInputGeom(end+1)=ftell(fid);
            end
        end

        if ~isempty(posInputGeom)
            geomslist = posInputGeom;
        else
            geomslist = posStandGeom;
        end
            
        if fl_createmovie || fl_extractall
            num=1;
            for pos=geomslist
                fseek(fid,pos,'bof');
                [ms0,status]=extrgeom(fid,1,0);
                if status
                    disp([worktitle,': ',lastwarn]);
                    fclose(fid);
                    continue
                end

                if fl_extractall
                    %---------
                    ms0.desc = [fnameshort '_#' num2str(num,'%02.0f')];

                    if ~fl_identmol
                        order=1:ms0.atomnum;
                    else
                        ms0 = createbondtable(ms0);
                        [ms0,status] = identmol(ms0,moltype);
                        if status
                          disp([worktitle,': ',lastwarn]);
                          fclose(fid);
                          continue
                        end
                        ms0 = calcproperties(ms0,moltype);

                        [ms0,desc]=checknewconf(ms0,workdb,desc);

                        order=[];
                        for I=1:numel(orderanchor)
                          order(end+1)=find(ms0.pind==orderanchor(I));
                        end

                        ms0=createzmt(ms0,order);
                        ms0=zmt2xyz(ms0,order);
                    end

                    if savemode.xyz
                        savemol(odir,ms0,0,order);
                    end
				    if savemode.zmt
				        savemol(odir,ms0,3,order);
					end
                    if savemode.gsxyz
                      savemolgs(odir,ms0,3,order,fullgtemplname); %Gaussian with XYZ
                    end
                    if savemode.gszmt
                      savemolgs(odir,ms0,5,order,fullgtemplname); %Gaussian with ZMT with varlist
                    end
                    %---------------    

                end
                if fl_createmovie || fl_createmovie
                    cla
                    plotmol(ms0,'b',0);
                    axis_saved=axis;
                    text((axis_saved(1)+axis_saved(2))/2,axis_saved(4)-0.05*(axis_saved(4)-axis_saved(3)),[fnameshort ': mode #' int2str(num)])
                    pause
                end
                num=num+1;
            end
        else %~ (fl_createmovie || fl_extractall)


            fseek(fid,geomslist(end),'bof');
            [ms0,status]=extrgeom(fid,1,0);
            if status
                disp([worktitle,': ',lastwarn]);
                fclose(fid);
                continue
            end

            %---------
            ms0.desc = fnameshort;

            if ~fl_identmol
                order=1:ms0.atomnum;
            else
                ms0 = createbondtable(ms0);
                [ms0,status] = identmol(ms0,moltype);
                if status
                  disp([worktitle,': ',lastwarn]);
                  fclose(fid);
                  continue
                end
                ms0 = calcproperties(ms0,moltype);

                [ms0,desc]=checknewconf(ms0,workdb,desc);

                order=[];
                for I=1:numel(orderanchor)
                  order(end+1)=find(ms0.pind==orderanchor(I));
                end

                ms0=createzmt(ms0,order);
                ms0=zmt2xyz(ms0,order);
            end

            if savemode.xyz
                savemol(odir,ms0,0,order);
            end
			if savemode.zmt
			    savemol(odir,ms0,3,order);
			end
            if savemode.gsxyz
                savemolgs(odir,ms0,3,order,fullgtemplname); %Gaussian with XYZ
            end
            if savemode.gszmt
                savemolgs(odir,ms0,5,order,fullgtemplname); %Gaussian with ZMT with varlist
            end
            %---------------    
       end
       fclose(fid);
%    catch
%       disp(['error: ',lasterr]);
%       fclose(fid);
%    end
end %for

diary off

