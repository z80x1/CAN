%rot15: analyze directory with xyz files for new conformations
%
% Version 1.0    
% Last modified  R O Zhurakivsky 2008-05-18
% Created        R O Zhurakivsky 2006-?-?

global pind

clear ms0 workdb desc filenames
format compact
atomsind


moltype=9
%usedpackage='xyz'
usedpackage='Gaussian'
theory='dftV2'
onlyoriginal=0  %#ok %process db with only original conformations
flwritefile=1   %#ok

%indirdef='D:\_diplom\work\0611-0706.dUrd_eng\070803_crystal\out';
indirdef=[CD.datadir filesep 'g.out\r9\r927\v4\xyz_TS' ]


if moltype==7
	  indir='r7.xyz';
      indir=[CD.xyzdir filesep indir]
      workdbname='r7';
elseif moltype==8
	  indir='r8z';
      indir=[CD.xyzdir filesep indir]
      workdbname='r8';
elseif moltype==11
	  indir='azacyt';
      indir=[CD.xyzdir filesep indir]
      workdbname='r11';
%elseif moltype==9
%	  indir='r9.fromNWChem';
%      indir=[CD.xyzdir filesep indir]
%      workdbname='r9';
else
    indir=indirdef
    workdbname=['r' int2str(moltype)]   %#ok
end


if usedpackage=='Gaussian'
  workdbname=[workdbname '_g'];
end
if ~strcmp(theory,'dft')
  workdbname=[workdbname '_' theory];
end
if onlyoriginal
    templ='_or';
    workdbname = [workdbname templ];
end

workdbname=[CD.dbdir filesep workdbname '.mat']

if exist(workdbname)==2
  load(workdbname,'workdb');
else
  workdb={};  
end


desc=char();
if ~isempty(workdb)
  decssize = numel(workdb(1).prop.sdesc);
  for i=1:numel(workdb)
      desc(i,1:decssize) = workdb(i).prop.sdesc(1:decssize);
%      sdescnew(i) = {workdb(i).prop.sdescnew};
%      worktitles(i) = {workdb(i).nwworktitle};
      filenames(i) = {lower(workdb(i).worktitle)};
  end
else
  filenames={};
end

badfilelist=cell(0);

sfiles = dir(strcat(indir,filesep,'*.xyz'));
numfiles = size(sfiles,1);

if numfiles==0
   error('No files to process found');
end

maxind=0;
blen=zeros(numfiles,1);
for i=1:numfiles

%    delete(gca);

    dlm=strfind(sfiles(i).name,'.');
    fnameshort = sfiles(i).name(1:(dlm-1));
    fnamefull = sfiles(i).name(1:(dlm(end)-1));

    if strcmpcellar(filenames,fnamefull)
      continue
    end

    ffname = strcat(indir,filesep,sfiles(i).name); 

    strcat('Loaded file: ',ffname)


	ms0=loadmolxyz(ffname,fnameshort,maxind);
   

    ms0 = createbondtable(ms0);
    [ms0,status] = identmol(ms0,moltype);
    if status
		disp(['error: ' lastwarn])
%        fclose(fid);
        continue
    end
    ms0 = createzmt(ms0);
  %    ms0=zmt2xyz(ms0);
    ms0 = calcproperties(ms0,moltype);

    ms0.GO=[];
    ms0.SP=[];
    ms0.DM=[];

    [ms0,desc]=checknewconf(ms0,workdb,desc);

    ms0.filename=fnamefull;
    ms0.worktitle=[];

    ms0.ZPE=[];
	ms0.gaussian.ZPE=[];
	ms0.freq.freq=0;
	ms0.comment={''};
	ms0.date=now;

if flwritefile 
      if isempty(workdb)
        workdb=ms0;
      else
        ms0=orderfields(ms0,workdb(1));
        workdb(end+1)=ms0;
      end  
      filenames(end+1) = {workdb(end).filename};
end
 
%    'Press any key.'
%    pause
end

save(workdbname,'workdb')

