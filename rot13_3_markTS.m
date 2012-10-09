%mark conformations with negative lowest vibrational mode freuency with new='T' mark (if not is done before)
%
% Version 1.0    
% Last modified  R O Zhurakivsky 2008-07-30
% Created        R O Zhurakivsky 2006-?-?

tic
clear 
format compact

global pind
atomsind

%----------------------------------------
moltype=412
usedpackage='Gaussian'
%theory='dftV2'
theory='dftV4'
%----------------------------------------


workdbname=['r' int2str(moltype)]
if usedpackage=='Gaussian'
  workdbname=[workdbname '_g'];
end
if ~strcmp(theory,'dft')
  workdbname=[workdbname '_' theory];
end
workdbname=[CD.dbdir filesep workdbname '.mat']

numchanged=0;

load(workdbname,'workdb')
recnum=numel(workdb);
for i=1:recnum
%disp(i)

    if isfield(workdb(i),'freq')

      if min(workdb(i).freq.freq)<0

        disp([workdb(i).prop.sdesc '(#' int2str(i) '): lowest freq = ' num2str(min(workdb(i).freq.freq)) ', new=' workdb(i).new])
        if workdb(i).new=='Y'
          workdb(i).new='T';
      disp(['Changed to TS: ' workdb(i).prop.sdesc '(#' int2str(i) ')'])
          numchanged = numchanged + 1;
        end
      end
    end
end

disp(['  Totally ' int2str(numchanged) ' conformations changed to TS.'])
%workdb = workdb_new;

if numchanged
  dlm=strfind(workdbname,'.');
  workdbnameold=[workdbname(1:dlm(end)-1) '~' workdbname(dlm(end):end)];
  copyfile(workdbname,workdbnameold);

  save(workdbname,'workdb')
else
  disp('Results are not saved!')
end

toc