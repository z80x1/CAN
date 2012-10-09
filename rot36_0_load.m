%rot36_AIMcorrel: analyze correlations for all Hbonds
%
% Version 1.0    
% Last modified  R O Zhurakivsky 2009-04-25
% Created        R O Zhurakivsky 2009-04-24

tic
clear 
format compact
atomsind


%---------------------------------
%test set
%moltype=[110 152 310 352 353 390 391 161 162] %#ok
%theory =[3   3   3   3   3   3   3   3   3];

%DNA & RNA nucleosides set
moltype=[9 13 15 12 16 140 240 340 540 ] %#ok
theory =[3  3  3  2  3   2   2   2   2 ];

%full set
%moltype=[9 151 14 13 15 17 19 21 270 521 11 12 16 140 240 340 540 110 152 310 352 353 390 391 161 162] %#ok
%theory =[3   2  3  3  3  3  3  3   2   3  2  2  3   2   2   2   2   3   3   3   3   3   3   3   3   3];

usedpackage='Gaussian' %#ok
 
theory=strcat ( 'dftV', int2str(theory'));
onlyoriginal=1;  % process db with only original conformations
%---------------------------------

%db.workdb=struct([]);
%db.moltype=zeros(1,numel(moltype));
%db.title=[];

for m_ind=1:numel(moltype)


    workdbname=['r' int2str(moltype(m_ind))] %#ok
    if strcmp(usedpackage,'Gaussian')
      workdbname=[workdbname '_g'];
    end
    if ~strcmp(theory,'dft')
      workdbname=[workdbname '_' theory(m_ind,:)];
    end
    if onlyoriginal
        templ='_or';
        workdbname = [workdbname templ];
    end
    workdbname=[CD.dbdir filesep workdbname '.mat'] %#ok

    load(workdbname,'workdb')
    recnum=numel(workdb);

    db(m_ind).workdb=workdb;
    db(m_ind).moltype=moltype(m_ind);
    db(m_ind).title=['r' int2str(moltype(m_ind)) '_' theory(m_ind,:)];

end
toc
