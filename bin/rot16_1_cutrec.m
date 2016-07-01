%cut conformations database to specified number of records
%
% Version 1.0    
% Last modified  R O Zhurakivsky 2008-05-19
% Created        R O Zhurakivsky 2006-?-?

workdbname='r12_g'  %#ok
workdbname=[CD.dbdir filesep workdbname '.mat'] %#ok

load(workdbname,'workdb')

A=workdb(1:148);
workdb=A;

%recnum=numel(workdb);
%for i=1:recnum
%  if isfield(workdb(i).gaussian,fieldname)
%    workdb(i).gaussian=rmfield(workdb(i).gaussian,fieldname);
%  end
%end
save(workdbname,'workdb')
