%rot29: sort XyyyZ conformations by bbb values
%
% Version 1.0    
% Last modified  R O Zhurakivsky 2008-05-19
% Created        R O Zhurakivsky 2006-?-?

tic
clear 
format compact

global pind
atomsind

moltype=8
usedpackage='Gaussian'
theory='mp2V2'
onlyoriginal=1;  % process db with only original conformations

workdbname=['r' int2str(moltype)]
if usedpackage=='Gaussian'
  workdbname=[workdbname '_g'];
end
if ~strcmp(theory,'dft')
  workdbname=[workdbname '_' theory];
end
if onlyoriginal
    templ='_or';
    workdbname = [workdbname templ];
end
workdbname=[CD.dbdir filesep workdbname '.mat']

load(workdbname,'workdb')

recnum=numel(workdb);

tailsconf=char(recnum,3);
sdesclen=numel(workdb(1).prop.sdesc);
sdesc={};
sdesc4sort=char(zeros(recnum,sdesclen));

for i=1:recnum

   tailsconf(i,1:3) = workdb(i).prop.sdesc(2:4);
   sdesc{i} = workdb(i).prop.sdesc;
   if sdesclen==5
     sdesc4sort(i,1:end) = [ workdb(i).prop.sdesc(5) workdb(i).prop.sdesc(1) workdb(i).prop.sdesc(2:4) ];
   else
     sdesc4sort(i,1:end) = workdb(i).prop.sdesc;
   end
end


tailsconfu=unique(tailsconf,'rows');

uniqtailnum=size(tailsconfu,1);
for i=1:uniqtailnum
  inds=find(sum(tailsconf==repmat(tailsconfu(i,:),recnum,1),2)==3);

  [XXX,II]=sortrows(sdesc4sort(inds,:));
  inds = inds(II);

  disp(sdesc(inds))
end

disp(['Total conformations: ' int2str(recnum)])
disp(['Total unique yyy combinations: ' int2str(uniqtailnum)])

toc
