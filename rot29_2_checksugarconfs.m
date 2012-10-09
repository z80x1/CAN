%rot29_2_checksugarconfs: compare ribose sugar confs with sugar confs in nucleosides
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
nucmoltype=[9 12 13 15 16]

usedpackage='Gaussian'
theory='dftV2'
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

sdesclen=numel(workdb(1).prop.sdesc);
sdescR={}; %ribose conformations
sdescR4sort=char(zeros(recnum,3));

for i=1:recnum
   sdescR{i} = workdb(i).prop.sdesc;
   sdescR4sort(i,1:end) = [ workdb(i).prop.sdesc(2:4) ];
end
clear workdb

nuclnum=numel(nucmoltype);
sdescNS=cell(40,nuclnum);  %nucleosides syn
sdescNS4sort=char(zeros(40,3,nuclnum)); 
sdescNA=cell(40,nuclnum);  %nucleosides anti
sdescNA4sort=char(zeros(40,3,nuclnum));
maxrecnum=0;
maxsyn=0; maxanti=0;
for jj=1:nuclnum

    moltype=nucmoltype(jj);
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
    maxrecnum=max(maxrecnum,recnum);

    sdesclen=numel(workdb(1).prop.sdesc);

    isyn=0;
    ianti=0;
    for i=1:recnum
        if workdb(i).prop.sdesc(end)=='A'
           ianti=ianti+1;
           sdescNA{ianti,jj} = workdb(i).prop.sdesc;
           sdescNA4sort(ianti,1:end,jj) = [ workdb(i).prop.sdesc(2:4) ];
        else
           isyn=isyn+1;
           sdescNS{isyn,jj} = workdb(i).prop.sdesc;
           sdescNS4sort(isyn,1:end,jj) = [ workdb(i).prop.sdesc(2:4) ];
        end
    end
    maxsyn=max(maxsyn,isyn);
    maxanti=max(maxanti,ianti);
    
    clear workdb
end
%sdescNS=sdescNS{1:maxsyn,jj};
%sdescNS4sort=sdescNS4sort(1:maxsyn,:,jj);
%sdescNA=sdescNA{1:maxanti,jj};
%sdescNA4sort=sdescNA4sort(1:maxanti,:,jj);

sdescNS=sdescNS(:);

toc
