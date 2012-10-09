%rot22_3: compare main torsion angles sections for selected molecules
%
% Version 1.0    
% Last modified  R O Zhurakivsky 2011-10-24
% Created        R O Zhurakivsky 2011-10-24

atomsind
format compact


%-------------------------------------------
%moltypes= [9 151 14 13 15 17 19 21 270 521 11 12 16 140 240 340 540 110 152 310 352 353 390 391 161 162] %#ok
%theories =[3   2  1  3  3  3  3  3   2   2  2  2  2   2   2   2   2   3   3   3   3   3   3   3   3   3];
moltypes= [15 16 161 162 390 391] %#ok
theories =[3 3 3 3 3 3];

usedpackage='Gaussian'  %#ok
%theory='dft'          %#ok
onlyoriginal=1  %#ok %process db with only original conformations

outfile = [CD.xlsdir filesep 'torsion_sect_comp' '.csv'];

%-------------------------------------------

theories=strcat ( 'dftV', int2str(theories'));

if numel(moltypes)>1
    fl_multiplemol=1;
else
    fl_multiplemol=0;
end

molnum = numel(moltypes);

prop = [];
% prop.beta = zeros(molnum,4);
% prop.gamma = zeros(molnum,4);
% prop.epsilon = zeros(molnum,4);
% prop.P = zeros(molnum,4);

T = {};
T(1,1) = {'molecule'};
T(1,2) = {'param'};
T(1,3) = {'beta'};
T(1,6) = {'gamma'};
T(1,9) = {'epsilon'};
T(1,12) = {'delta'};
T(1,14) = {'P'};
T(1,17) = {'chi'};
T(2,3:20) = [{'g+'} {'t'} {'g-'} {'g+'} {'t'} {'g-'} {'g+'} {'t'} {'g-'} {'g+'} {'t'} ...
             {'N'} {'E'} {'S'} {'syn'} {'high-syn'} {'anti'} {'high-anti'}];

row=3;

for mind=1:molnum
    moltype=moltypes(mind);
    theory=theories(mind,:);

    if any(moltype==[8])
        ringtypes=['b'; 'g'; 'd'; 'e'; 'p'];
    elseif any(moltype==[7])
        ringtypes=['t'; 'n'; 'b'; 'g'; 'd'; 'e'; 'p'];
    elseif any(moltype==[112,212,312,412,512])
        ringtypes=['b'; 'g'; 'c'; 'p';];
    elseif any(moltype==[11,15,16,9,12,13,521,161,162,110,152,310,352,353,390,391,413,523,524])
        ringtypes=['b'; 'g'; 'd'; 'e'; 'c'; 'p'];
    %    ringtypes=['b'; 'g'; 'd'; 'e'; 'c'; 'p'; '2'; '3'; '4'];
    elseif any(moltype==[140,240,340,540,14])
        ringtypes=['t'; 'n'; 'b'; 'g'; 'd'; 'e'; 'c'; 'p'];
    end

    %inclset={'AabcA'; 'EabcA'; 'JabcA'}

    workdbname=['r' int2str(moltype)];
    if strcmp(usedpackage,'Gaussian')
      workdbname=[workdbname '_g'];
    end
    if ~strcmp(theory,'dft')
      workdbname=[workdbname '_' theory];
    end
    if onlyoriginal
        workdbname = [workdbname '_or'];
    end
    workdbfname=[CD.dbdir filesep workdbname '.mat'] %#ok

    if exist(workdbfname,'file')==2
      load(workdbfname,'workdb');
%      disp(['Loaded ' workdbname]);
    else
      error('rot22_2:dbnexist','Database doesn''t exists!');
    end

    %ringtype='p' % b, g, d, e, c, p, t, n
    %clear tbeta tgamma tepsilon tdelta tchi tteta teta  Pdeg
    [tbeta,tgamma,tdelta,tepsilon,tchi,tteta,teta,Pdeg]=deal(zeros(0));
    
    recnum=numel(workdb);
    for i=1:recnum
      if workdb(i).new=='Y'

        ms0=workdb(i);

        if exist('inclset','var') && ~isempty(inclset)
           if isempty(intersect(ms0.prop.sdesc,inclset))
             continue
           end
        end

        tbeta(end+1)=ms0.prop.tbeta;
        tgamma(end+1)=ms0.prop.tgamma;

        if isfield(ms0.prop, 'tdelta')
	        tdelta(end+1)=ms0.prop.tdelta;
	    end
        if isfield(ms0.prop, 'tepsilon')
	        tepsilon(end+1)=ms0.prop.tepsilon;
	    end
        if isfield(ms0.prop,'tchi')
            tchi(end+1)=ms0.prop.tchi;
        end
        if isfield(ms0.prop,'tteta')
            tteta(end+1)=ms0.prop.tteta;
        end
        if isfield(ms0.prop,'teta')
            teta(end+1)=ms0.prop.teta;
        end
    %    P(end+1)=ms0.prop.P;
        Pdeg(end+1)=ms0.prop.Pdeg;

      end
    end
    
    T(row,1) = {['r' int2str(moltype)]};    
    T(row,2) = {'min'};    
    T(row+1,2) = {'max'};    
    T(row+2,2) = {'mean'};    
    T(row+3,2) = {'std'};    

    col=3;
    a = TorsionStat(tbeta(tbeta>0 & tbeta<120));
    T(row:row+3,col) = a.stat;
    a = TorsionStat([tbeta(tbeta>120) tbeta(tbeta<-120)+360]);
    T(row:row+3,col+1) = a.stat;
    a = TorsionStat(tbeta(tbeta>-120 & tbeta<0));
    T(row:row+3,col+2) = a.stat;
    col = col+3;

    a = TorsionStat(tgamma(tgamma>0 & tgamma<120));
    T(row:row+3,col) = a.stat;
    a = TorsionStat([tgamma(tgamma>120) tgamma(tgamma<-120)+360]);
    T(row:row+3,col+1) = a.stat;
    a = TorsionStat(tgamma(tgamma>-120 & tgamma<0));
    T(row:row+3,col+2) = a.stat;
    col = col+3;

    a = TorsionStat(tepsilon(tepsilon>0 & tepsilon<120));
    T(row:row+3,col) = a.stat;
    a = TorsionStat([tepsilon(tepsilon>120) tepsilon(tepsilon<-120)+360]);
    T(row:row+3,col+1) = a.stat;
    a = TorsionStat(tepsilon(tepsilon>-120 & tepsilon<0));
    T(row:row+3,col+2) = a.stat;
    col = col+3;
    
    a = TorsionStat(tdelta(tdelta>0 & tdelta<120));
    T(row:row+3,col) = a.stat;
    a = TorsionStat([tdelta(tdelta>120) tdelta(tdelta<-120)+360]);
    T(row:row+3,col+1) = a.stat;
    col = col+2;
    
    a = TorsionStat([Pdeg(Pdeg>=315)-360 Pdeg(Pdeg<=45)]); %N
    T(row:row+3,col) = a.stat;
    a = TorsionStat(Pdeg(Pdeg>45 & Pdeg<135)); %E
    T(row:row+3,col+1) = a.stat;
    a = TorsionStat(Pdeg(Pdeg>=135 & Pdeg<=225)); %S
    T(row:row+3,col+2) = a.stat;
%    a = TorsionStat(Pdeg(Pdeg>225 & Pdeg<315)); %W
%    T(row:row+3,col+3) = a.stat;
    col = col+3;
    
        %as in Saenger (p.33-35)    60deg<chi<=120deg - high-syn (V),  -60deg<chi<=60deg - syn (S)
        %               -150deg<chi<=-60deg - high-anti (B),  -180deg<chi<=-150deg,120deg<chi<=180deg - anti (B)
    a = TorsionStat(tchi(tchi>-60 & tchi<=60)); %syn
    T(row:row+3,col) = a.stat;
    a = TorsionStat(tchi(tchi>60 & tchi<=120)); %high-syn
    T(row:row+3,col+1) = a.stat;
    a = TorsionStat([tchi(tchi<-150)+360 tchi(tchi>120)]); %anti
    T(row:row+3,col+2) = a.stat;
    a = TorsionStat(tchi(tchi>-150 & tchi<=-60)); %high-anti
    T(row:row+3,col+3) = a.stat;

    
    row = row+4;
end %for mind

cell2csv(outfile,T,';',2003);       
disp('Done!')


