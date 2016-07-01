%build new database with only 'Y' conformations
%
% Version 1.0    
% Last modified  R O Zhurakivsky 2008-08-01
% Created        R O Zhurakivsky 2006-?-?

tic
atomsind

%----------------------------------------
moltype=16
theory='dftV3'  %#ok
%theory='mp2V2'  %#ok
usedpackage='Gaussian'  %#ok
%----------------------------------------

workdbname=[CD.dbdir filesep 'r' int2str(moltype)];
if strcmp(usedpackage,'Gaussian')
  workdbname=[workdbname '_g'];
end
if ~strcmp(theory,'dft')
  workdbname=[workdbname '_' theory];
end
newworkdbname=[workdbname '_or.mat'] %#ok
workdbname=[workdbname '.mat'] %#ok

load(workdbname,'workdb')

recnum=numel(workdb);

tmpdb={};
for i=1:recnum
  if workdb(i).new=='Y'

    if isempty(tmpdb)
        tmpdb=workdb(i);
    else
        tmpdb(end+1)=workdb(i);
    end
  end
end


workdb=tmpdb;
dlm=strfind(newworkdbname,'.');
newworkdbnamebkp=[newworkdbname(1:dlm(end)-1) '~' newworkdbname(dlm(end):end)];
if exist(newworkdbname,'file')
    copyfile(newworkdbname,newworkdbnamebkp);
end
save(newworkdbname,'workdb')

clear tmpdb

toc
