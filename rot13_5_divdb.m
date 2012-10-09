%rot13_5: divide db to two parts: main data and FC
%
% Version 1.0    
% Last modified  R O Zhurakivsky 2010-07-29
% Created        R O Zhurakivsky 2010-07-29

format compact
%!!!clear 
atomsind

%----------------------------
moltype=950  %#o
theory='dftV3'  %#ok
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


tic
load(workdbname,'workdb')

recnum=numel(workdb);
sdesc={};

clear workdb_main;
clear workdb_fc;
clear workdb_props;

fields = fieldnames(workdb);
for fi=1:numel(fields)
        fname=fields{fi};
        
        if strcmp(fname,'prop')
            for i=1:recnum
                ffields = fieldnames(workdb(i).prop);
            	for ffi=1:numel(ffields)
                    fname=ffields{ffi};
                    workdb_props(i).(fname)=workdb(i).prop.(fname)    ;
                end
            end

        elseif strcmp(fname,'freq')
            for i=1:recnum
                ffields = fieldnames(workdb(i).freq);
            	for ffi=1:numel(ffields)
                    fname=ffields{ffi};
                    workdb_fc(i).(fname)=workdb(i).freq.(fname)    ;
                end
            end

        else
            for i=1:recnum
                workdb_main(i).(fname)=workdb(i).(fname)    ;
            end
        end

end

workdb = workdb_main;
clear workdb_main;

dlm=strfind(workdbname,'.');
workdbnameold=[workdbname(1:dlm(end)-1) '~' workdbname(dlm(end):end)];
workdbname_fc=[workdbname(1:dlm(end)-1) '_fc' workdbname(dlm(end):end)];
workdbname_props=[workdbname(1:dlm(end)-1) '_props' workdbname(dlm(end):end)];

movefile(workdbname,workdbnameold);
save(workdbname,'workdb')
save(workdbname_fc,'workdb_fc')
clear workdb_fc;
save(workdbname_props,'workdb_props')
clear workdb_props;


toc
