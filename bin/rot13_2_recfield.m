%recalculate some DB field
%
% Version 1.0    
% Last modified  R O Zhurakivsky 2008-05-19
% Created        R O Zhurakivsky 2006-?-?

tic
clear 
format compact

atomsind

moltype=8   %#ok
usedpackage='Gaussian'  %#ok
theory='dftV2'  %#ok


workdbname=['r' int2str(moltype)]   %#ok
if strcmp(usedpackage,'Gaussian')
  workdbname=[workdbname '_g'];
end
if ~strcmp(theory,'dft')
  workdbname=[workdbname '_' theory];
end
workdbname=[CD.dbdir filesep workdbname '.mat'] %#ok


load(workdbname,'workdb')
recnum=numel(workdb);
for i=1:recnum
%disp(i)

    fnameshort = workdb(i).desc;
    dlm2=strfind(fnameshort,'_');
    if ~isempty(dlm2)
	if dlm2==numel(fnameshort)
	    fnameshort=fnameshort(1:end-1)  %#ok
	else	
            fnameshort=fnameshort(dlm2(end)+1:end)  %#ok
	end
	workdb(i).desc = fnameshort;
    end


end

%workdb = workdb_new;

dlm=strfind(workdbname,'.');
workdbnameold=[workdbname(1:dlm(end)-1) '~' workdbname(dlm(end):end)];
copyfile(workdbname,workdbnameold);

save(workdbname,'workdb')

toc

