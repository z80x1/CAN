%initialize some record field in conformations database
%
% Version 1.0    
% Last modified  R O Zhurakivsky 2008-05-19
% Created        R O Zhurakivsky 2006-?-?

atomsind

workdbname='r11_g_dftV3';
workdbfname=[CD.dbdir filesep workdbname '.mat']    %#ok

fieldname='MP2_631G';

load(workdbfname,'workdb')

recnum=numel(workdb);
for i=1:recnum
  if isfield(workdb(i).gaussian,fieldname)
    workdb(i).gaussian=rmfield(workdb(i).gaussian,fieldname);
  end
end
%workdb=workdbnew;

movefile(workdbfname,[CD.dbdir filesep workdbname '~.mat']);
save(workdbfname,'workdb')
