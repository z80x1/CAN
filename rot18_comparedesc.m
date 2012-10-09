%compare descriptions of two sets of conformations
%
% Version 1.0    
% Last modified  R O Zhurakivsky 2008-05-19
% Created        R O Zhurakivsky 2006-?-?

if 1
  workdbname=[CD.dbdir filesep 'r9_g.mat']  %#ok
  load(workdbname,'workdb')

  recnum=numel(workdb);
  set1desc=cell(1,recnum);
  for i=1:recnum
  %    set1desc(end+1)={workdb(i).filename}; %names of result files I have
      set1desc(i)={workdb(i).prop.sdesc};
      new1(i)=workdb(i).new;
  end
  [set1desc,ind1,xxx]=unique(set1desc);
  new1=new1(ind1);
  notnewconf1=set1desc(new1=='T')   %#ok

  workdbname=[CD.dbdir filesep 'r13_g.mat'] %#ok
  load(workdbname,'workdb')

  recnum=numel(workdb);
  set2desc=cell(1,recnum);
  for i=1:recnum
      set2desc(i)={workdb(i).prop.sdesc};
      new2(i)=workdb(i).new;
  end
  [set2desc,ind2,xxx]=unique(set2desc);
  new2=new2(ind2);
  notnewconf2=set2desc(new2=='T')   %#ok

  numels1=numel(set1desc)   %#ok
  numels2=numel(set2desc)   %#ok

  res12=setdiff(set2desc,set1desc)'  %#ok
  numel(res12)

  res21=setdiff(set1desc,set2desc)'  %#ok
  numel(res21)
end

if 0

  indir=[CD.xyzdir filesep 'rd01'];
  sfiles = dir(strcat(indir,filesep,'*.gjf'));
  numfiles = size(sfiles,1);
  for i=1:numfiles
    dlm1=strfind(sfiles(i).name,'_');
    if isempty(dlm1) 
  	dlm1=0;
    end	
    dlm2=strfind(sfiles(i).name,'.');
    set3desc(i) = {sfiles(i).name((dlm1(end)+1):(dlm2-1))}; %complete set of names
  end
  %set3desc

  indir=[CD.xyzdir filesep 'rd02'];
  sfiles = dir(strcat(indir,filesep,'*.gjf'));
  numfiles = size(sfiles,1);
  for i=1:numfiles
    dlm1=strfind(sfiles(i).name,'_');
    if isempty(dlm1) 
  	dlm1=0;
    end	
    dlm2=strfind(sfiles(i).name,'.');
    set4desc(i) = {sfiles(i).name((dlm1(end)+1):(dlm2-1))}; %complete set of names
  end
  %set4desc

  indir=[CD.xyzdir filesep 'rd03'];
  sfiles = dir(strcat(indir,filesep,'*.gjf'));
  numfiles = size(sfiles,1);
  for i=1:numfiles
    dlm1=strfind(sfiles(i).name,'_');
    if isempty(dlm1) 
  	dlm1=0;
    end	
    dlm2=strfind(sfiles(i).name,'.');
    set5desc(i) = {sfiles(i).name((dlm1(end)+1):(dlm2-1))}; %complete set of names
  end
  %set5desc

  numels3=numel(set3desc)   %#ok
  numels4=numel(set4desc)   %#ok
  numels5=numel(set5desc)   %#ok

  res=setdiff(set3desc,union(set4desc,set5desc))'  %#ok
  numel(res)

  indir=[CD.xyzdir filesep 'rd01'];
  for ii=1:numel(res)
    fname=[ indir filesep '*' res{ii} '*'];
    movefile(fname,[indir '_2'])
  end

end
