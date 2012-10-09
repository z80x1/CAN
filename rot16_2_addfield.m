%initialize some record field in conformations database
%
% Version 1.0    
% Last modified  R O Zhurakivsky 2008-05-19
% Created        R O Zhurakivsky 2006-?-?

atomsind

workdbname='r12_g_dft420_or'    %#ok
workdbname=[CD.dbdir filesep workdbname '.mat'] %#ok


load(workdbname,'workdb')

recnum=numel(workdb);
for i=1:recnum
    workdb(i).masses=zeros(size(workdb(i).pind));
end
%workdb=workdbnew;
save(workdbname,'workdb')
