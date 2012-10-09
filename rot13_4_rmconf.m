%rot13_4: remove conformations from db from specified list
%
% Version 1.0    
% Last modified  R O Zhurakivsky 2008-05-19
% Created        R O Zhurakivsky 2006-?-?

format compact
clear 
atomsind

%----------------------------
moltype=240  %#ok
theory='dftV2'  %#ok
onlyoriginal=0;  % process db with only original conformations

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
workdbname=[workdbname '.mat']  %#ok


workname='r24013' %#ok
fileslist.name=[CD.xyzdir filesep workname filesep 'fileslist'];
fileslist  %#ok

tic
load(workdbname,'workdb')

recnum=numel(workdb);
sdesc={};
worktitle={};
filename={};
for i=1:recnum

    sdesc(i) = {workdb(i).prop.sdesc};
    worktitle(i) = {workdb(i).worktitle};
    filename(i) = {workdb(i).filename};
end

rmfileset = textread(fileslist.name,'%s');
[rem,I]=setdiff(filename,rmfileset );
workdb_new=workdb(sort(I));

workdb = workdb_new;

dlm=strfind(workdbname,'.');
workdbnameold=[workdbname(1:dlm(end)-1) '~' workdbname(dlm(end):end)];
copyfile(workdbname,workdbnameold);

save(workdbname,'workdb')

toc
