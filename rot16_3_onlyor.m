%build new database with only 'Y' conformations
%
% Version 1.0    
% Last modified  R O Zhurakivsky 2008-08-01
% Created        R O Zhurakivsky 2006-?-?

% 2010-0820 it is nonoptimal to load all three parts of DB at one time, need to be rewritten

tic
atomsind

%----------------------------------------
moltype=8
theory='dftV3bis'  %#ok
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
orworkdbname=[workdbname '_or.mat'] %#ok
workdbname=[workdbname '.mat'] %#ok

load(workdbname,'workdb')

  fl_dbprops_exist=0;
  %zhr100820
  dlm=strfind(workdbname,'.');
  workdbname_props=[workdbname(1:dlm(end)-1) '_props' workdbname(dlm(end):end)];
  if exist(workdbname_props,'file')==2
    clear workdb_props;
    load(workdbname_props,'workdb_props');
    fl_dbprops_exist=1;
  end
  fl_dbfc_exist=0;
  dlm=strfind(workdbname,'.');
  workdbname_fc=[workdbname(1:dlm(end)-1) '_fc' workdbname(dlm(end):end)];
  if exist(workdbname_fc,'file')==2
    clear workdb_fc;
    load(workdbname_fc,'workdb_fc');
    fl_dbfc_exist=1;
  end

recnum=numel(workdb);

tmpdb={};
tmpdb_props={};
tmpdb_fc={};
for i=1:recnum
  if workdb(i).new=='Y'

    if isempty(tmpdb)
        tmpdb=workdb(i);
        if (fl_dbprops_exist)
		  tmpdb_props=workdb_props(i);
        end
        if (fl_dbfc_exist)
		  tmpdb_fc=workdb_fc(i);
        end
    else
        tmpdb(end+1)=workdb(i);
        if (fl_dbprops_exist)
		  tmpdb_props(end+1)=workdb_props(i);
        end
        if (fl_dbfc_exist)
		  tmpdb_fc(end+1)=workdb_fc(i);
        end
    end
  end
end


workdb=tmpdb;
clear tmpdb;
dlm=strfind(orworkdbname,'.');
orworkdbnamebkp=[orworkdbname(1:dlm(end)-1) '~' orworkdbname(dlm(end):end)];
orworkdbname_fc=[orworkdbname(1:dlm(end)-1) '_fc' orworkdbname(dlm(end):end)];
if exist(orworkdbname,'file')
    movefile(orworkdbname,orworkdbnamebkp);
end
save(orworkdbname,'workdb');
clear workdb;

if (fl_dbprops_exist)
  workdb_props=tmpdb_props;
  clear tmpdb;
  orworkdbname_props=[orworkdbname(1:dlm(end)-1) '_props' orworkdbname(dlm(end):end)];
  orworkdbnamebkp_props=[orworkdbname(1:dlm(end)-1) '_props~' orworkdbname(dlm(end):end)];
  if exist(orworkdbname_props,'file')
      movefile(orworkdbname_props,orworkdbnamebkp_props);
  end
  save(orworkdbname_props,'workdb_props')
end
clear workdb_props;

if (fl_dbfc_exist)
  workdb_fc=tmpdb_fc;
  clear tmpdb;
  orworkdbname_fc=[orworkdbname(1:dlm(end)-1) '_fc' orworkdbname(dlm(end):end)];
  orworkdbnamebkp_fc=[orworkdbname(1:dlm(end)-1) '_fc~' orworkdbname(dlm(end):end)];
  if exist(orworkdbname_fc,'file')
      movefile(orworkdbname_fc,orworkdbnamebkp_fc);
  end
  save(orworkdbname_fc,'workdb_fc')
end
clear workdb_fc;

toc
