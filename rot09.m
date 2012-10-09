%create excel file with molecule properties
%
% Version 1.0    
% Last modified  R O Zhurakivsky 2008-12-15
% Created        R O Zhurakivsky 2005-09-?

%clear 
format compact
time0 = cputime;
tic
pindsdef
atomsind


%--------- data section -------------------------------
T=298.15         %#ok Temperature
moltype = 7  %#ok
%usedpackage='NWChem'
%usedpackage='xyz'
usedpackage='Gaussian'  %#ok
%theory='mp2ccpVDZ'
theory='dftV3bis'  %#ok
%theory='dftV2'  %#ok
%theory='pair1003'  %#ok
%theory='cryst'  %#ok
onlyoriginal=1  %#ok %process db with only original conformations
dbsuffix=''
flwritefile=1   %#ok

version='01'    %#ok

indir = CD.dbdir;
%indir = 'E:\дезокситимидин\гауссиан\Newfolder';

outdir = CD.xlsdir;
%outdir = indir;
%--------- data section end ---------------------------

xlsfile = ['r' int2str(moltype)];
workdbname=['r' int2str(moltype)]   %#ok

if strcmp(usedpackage,'Gaussian')
  xlsfile = [xlsfile '_g'];
  workdbname=[workdbname '_g'];
end
if ~strcmp(theory,'dft')
  workdbname=[workdbname '_' theory];
  xlsfile=[xlsfile '_' theory];
end
if uint16(T)~=298
  xlsfile = [xlsfile '_' int2str(T)];
end

if ~isempty(version)
  xlsfile=[xlsfile '_' version];
end

if onlyoriginal
    templ='_or';

    xlsfile = [xlsfile templ];
    workdbname = [workdbname templ];

end
if ~isempty(dbsuffix)
    xlsfile = [xlsfile '_' dbsuffix];
    workdbname = [workdbname '_' dbsuffix];
end

xlsfile=[outdir filesep xlsfile '.xls']  %#ok
workdbname=[indir filesep workdbname '.mat'] %#ok

if exist(workdbname,'file')==2
  load(workdbname,'workdb');
else
  error(['File doesn''t exist: ' workdbname])
end
recnum=numel(workdb);

if recnum==0
    disp('Exited. Database is empty.');
    return
end    
if flwritefile
  disp('Have you sorted XLS file in record accending order? ');
  answer=input('Are you sure? [1]/0');
  if ~answer 
    return
  end
end


[ZPE,TEC,EEC,GEC]=deal([]);
itableHbond=[0 0];

Hbonddatadesc={};
if isfield(workdb(1).prop,'Hbonddatadesc')
  for i=1:recnum
      Hbonddatadesc=union(Hbonddatadesc,workdb(i).prop.Hbonddatadesc);
  end
end
Hbonddata=zeros(recnum,numel(Hbonddatadesc));
freq.freq=zeros(recnum,workdb(1).atomnum*3-6);
freq.intensity=zeros(recnum,workdb(1).atomnum*3-6);
freq.forceconst=zeros(recnum,workdb(1).atomnum*3-6);
%[HbondLMind,HbondVMind]=deal(zeros(recnum,numel(workdb(1).prop.HbondHind)));

if strcmp(usedpackage,'Gaussian')
  fdescgausenergy = {};
  for i=1:recnum
    fdescgausenergyadd=fieldnames(workdb(i).gaussian);
    fdescgausenergy=unique([fdescgausenergy;fdescgausenergyadd]);
  end
  fdescgausenergy_all=setdiff(fdescgausenergy',[{'stpoint'} {'DM'} {'ZPE'} {'TEC'},{'EEC'},{'GEC'},{'DMv'},{'imoments'},{'mcharge'},{'T'}]);

  fdescgausenergyGECcorr={};
  for i=1:numel(fdescgausenergy_all)
    buf=fdescgausenergy_all(i);
    fdescgausenergyGECcorr(i)={[buf{1} ' (GEC corrected)']};
  end
end
%find most accurate energy : MP2 is preferable then field with longest name is chosen
if ~isempty(fdescgausenergy_all)
    inds=strmatch('MP2',fdescgausenergy_all);
    if isempty(inds) 
      inds=numel(fdescgausenergy_all);
    end
    fnsize = repmat(0,1,numel(inds));
    for ii=1:numel(inds)
      fnsize(ii)=numel(fdescgausenergy_all{inds(ii)});
    end
    [XXX,I]=max(fnsize);
    mostaccenergyfn=fdescgausenergy_all{inds(I)};
    mostaccenergyfnind=strcmpcellar(fdescgausenergy_all,mostaccenergyfn);
end

%--distance & angles & torsion matrixes
%[pindsorted,pindsortii]=sort(pind.labels);
[distmatrorder,distmatrordii] = sortrows(workdb(1).prop.distmatrpinds);
distmatrdesc=cell(size(distmatrordii));
ind=0;
for i=distmatrordii'
  ind=ind+1;
  distmatrdesc(ind)={[pind.labels{workdb(1).prop.distmatrpinds(i,1)}, pind.labels{workdb(1).prop.distmatrpinds(i,2)}]};
end
[anglematrorder,anglematrordii] = sortrows(workdb(1).prop.anglematrpinds);
anglematrdesc=cell(size(anglematrordii));
ind=0;
for i=anglematrordii'
  ind=ind+1;
  anglematrdesc(ind)={[pind.labels{workdb(1).prop.anglematrpinds(i,1)} pind.labels{workdb(1).prop.anglematrpinds(i,2)}...
               pind.labels{workdb(1).prop.anglematrpinds(i,3)}]};
end
[tormatrorder,tormatrordii] = sortrows(workdb(1).prop.tormatrpinds);
tormatrdesc=cell(size(tormatrordii));
ind=0;
for i=tormatrordii'
  ind=ind+1;
  tormatrdesc(ind)={[pind.labels{workdb(1).prop.tormatrpinds(i,1)} pind.labels{workdb(1).prop.tormatrpinds(i,2)} ...
             pind.labels{workdb(1).prop.tormatrpinds(i,3)} pind.labels{workdb(1).prop.tormatrpinds(i,4)}]};
end
%--end

[hta0,hta1,hta2,hta3,hta4,hta5,tbeta,tgamma,tdelta,tepsilon,tchi,teta,tteta,tkappa1,tkappa2,...
 rglyc,tau01,tau02,tau03,tau04,tau05,tau06,tau10,tau11,tau12,maxtorinbase,...
 tnu0,tnu1,tnu2,tnu3,tnu4,rnu0,rnu1,rnu2,rnu3,rnu4,dplC2,dplC3,...
 aO4C1C2,aC1C2C3,aC2C3C4,aC3C4O4,aC4O4C1,sumofnu,Pdeg,numax,O5Hn,rO5H5,rO3H3,...
 SPenergy,GOenergy,nbonds,DM,GOsteps ] = ...
        deal(repmat(0,1,recnum));
torbpl = deal(repmat(0,recnum,1));
[desc,confdesc,sdesc,sdescnew,filename,worktitle,new,ldate,comment,torbpldesc,HbondHind] = ...
    deal(repmat({},1,recnum));
gausenergy = repmat(0,numel(fdescgausenergy_all),recnum);
[Hind,Hbondpindlabels,HbondLMind,HbondVMind,HbondVMintensity,HbondLMfreq,HbondVMfreq] = ...
    deal(repmat(0,recnum,numel(strcmpcellar(workdb(1).labels,'H'))));
HbondFC=[];
Hbondindprint=[];
Hbondpindlabels = repmat({},numel(strcmpcellar(workdb(1).labels,'H')),recnum);
distmatr = repmat(0,recnum,numel(workdb(1).prop.distmatr));
anglematr = repmat(0,recnum,numel(workdb(1).prop.anglematr));
tormatr = repmat(0,recnum,numel(workdb(1).prop.tormatr));

for i=1:recnum

    ms0=workdb(i);
    
    desc(i) = {ms0.desc};
    confdesc(i) = ms0.prop.confdesc;
    sdesc(i) = {ms0.prop.sdesc};
    sdescnew(i) = {ms0.prop.sdescnew};
    if strcmp(usedpackage,'Gaussian')
      if ms0.filename(1)=='='
          buf = ['''' ms0.filename];
      else
          buf = ms0.filename;
      end
      filename(i) = {buf};
      if ms0.worktitle 
          if ms0.worktitle(1)=='='
              buf = ['''' ms0.worktitle];
          else
              buf = ms0.worktitle;
          end
      else
          buf = '';
      end
      worktitle(i) = {buf};
    else
      filename(i) = {ms0.nwfilename};
      worktitle(i) = {ms0.nwworktitle};
    end
    new(i) = {ms0.new};

    if ~isfield(ms0,'date') || isempty(ms0.date)
        ldate(i) = {''};
    else
        ldate(i) = {datestr(ms0.date)};
    end
    for j=1:numel(ms0.comment)
      commnext=ms0.comment(j);
      if j==1
        comm = commnext{1};
      else
        comm = [comm '; ' commnext{1}];
      end   
    end
    comment(i)={comm};

    hta2(i) = ms0.prop.hta2;
    hta3(i) = ms0.prop.hta3;
    if isfield(ms0.prop,'hta4')
        hta4(i) = ms0.prop.hta4;
    else
        hta4(i) = NaN;
    end

    if any(moltype==[7,14]) || mod(moltype,100)==40 %non 2'deoxy molecules
        hta5(i) = ms0.prop.hta5;
    end    

    hta0(i) = ms0.prop.hta0;
    
    tbeta(i) = ms0.prop.tbeta;
    tgamma(i) = ms0.prop.tgamma;

    if isfield(ms0.prop,'tdelta');
        tdelta(i) = ms0.prop.tdelta;
    else
        tdelta(i) = NaN;
    end
    if isfield(ms0.prop,'tepsilon');
        tepsilon(i) = ms0.prop.tepsilon;
    else
        tepsilon(i) = NaN;
    end

    if any(moltype==[7,14]) || mod(moltype,100)==40
        teta(i) = ms0.prop.teta;
        tteta(i) = ms0.prop.tteta;
        tkappa1(i) = ms0.prop.tkappa1;
    end
    if any(moltype==[11,9,12,13,14,17,18,19,21]) || any(floor(moltype/100)==[2,4,5]) %Cyd,Thd,Urd
        tchi(i) = ms0.prop.tchi;
        rglyc(i) = ms0.prop.rC1N1;
    elseif any(moltype==[15,16]) || any(floor(moltype/100)==[1,3]) %Ade,Guo
        tchi(i) = ms0.prop.tchi;
        rglyc(i) = ms0.prop.rC1N9;
    end

    if ~any(moltype==[7,8])
        tau01(i) = ms0.prop.tau01; %angle between C1'N1 and aminobase plane (previously aglycoutplane)
        maxtorinbase(i) = ms0.prop.maxtorinbase;


        torbpl1 = ms0.prop.torbpl;
        torbpl2=[];
        
        if i==1
        torbplind = ms0.prop.torbplind;
        jend = numel(torbplind);
        for j=1:jend
            j1=j-1; if j1==0, j1=jend; end;
            j3=j+1; if j3>jend, j3=j3-jend; end;
            j4=j+2; if j4>jend, j4=j4-jend; end;
            s1=pind.labels{ms0.pind(torbplind(j1))};
            s2=pind.labels{ms0.pind(torbplind(j))};
            s3=pind.labels{ms0.pind(torbplind(j3))};
            s4=pind.labels{ms0.pind(torbplind(j4))};
            torbpldesc(j) = {[s1(2:end) s2(2:end) s3(2:end) s4(2:end)]};
        end
        end

        if any(moltype==[15,16]) || any(floor(moltype/100)==[1,3]) %Ade,Guo

            torbpl2 = ms0.prop.torbpl2;

            if i==1
            torbpl2ind = ms0.prop.torbpl2ind;
            jend2 = numel(torbpl2ind);
            for j=1:jend2
                j1=j-1; if j1==0, j1=jend2; end;
                j3=j+1; if j3>jend2, j3=j3-jend2; end;
                j4=j+2; if j4>jend2, j4=j4-jend2; end;
                s1=pind.labels{ms0.pind(torbpl2ind(j1))};
                s2=pind.labels{ms0.pind(torbpl2ind(j))};
                s3=pind.labels{ms0.pind(torbpl2ind(j3))};
                s4=pind.labels{ms0.pind(torbpl2ind(j4))};

                torbpldesc(j+jend) = {[s1(2:end) s2(2:end) s3(2:end) s4(2:end)]};
            end
            end
        end

        torbplsum = [torbpl1 torbpl2];
        inum = numel(torbplsum);
        torbpl(i,1:inum) = torbplsum;

        tau04(i) = ms0.prop.tau04;
        tau05(i) = ms0.prop.tau05;
        tau06(i) = ms0.prop.tau06;

    end
        
    if isfield(ms0.prop,'tau10')
        tau10(i) = ms0.prop.tau10; %angle between plane of NH2 and nucleicbase C2C4C5 plane 
    else tau10(i)=NaN; 
    end
    if isfield(ms0.prop,'tau11')
        tau11(i) = ms0.prop.tau11; %out of aminobase plane angle of H41 atom 
    else tau11(i)=NaN; 
    end
    if isfield(ms0.prop,'tau12')
        tau12(i) = ms0.prop.tau12; %out of aminobase plane angle of H42 atom 
    else tau12(i)=NaN; 
    end
    if isfield(ms0.prop,'tau02')
        tau02(i) = ms0.prop.tau02; %angle between C4N (or C2N for Guo) and aminobase least square fit plane
    else tau02(i)=NaN; 
    end
    if isfield(ms0.prop,'tau03')
        tau03(i) = ms0.prop.tau03; %angle between C4N (or C2N for Guo) and aminogroup plane NH2
    else tau03(i)=NaN; 
    end

    dplC2(i) = ms0.prop.dplC2;
    dplC3(i) = ms0.prop.dplC3;

    if isfield(ms0.prop,'tkappa2')
        tkappa2(i) = ms0.prop.tkappa2;
    else
        tkappa2(i) = NaN;
    end
    
    tnu0(i) = ms0.prop.tnu0;
    tnu1(i) = ms0.prop.tnu1;
    tnu2(i) = ms0.prop.tnu2;
    tnu3(i) = ms0.prop.tnu3;
    tnu4(i) = ms0.prop.tnu4;
    rnu0(i) = ms0.prop.rnu0;
    rnu1(i) = ms0.prop.rnu1;
    rnu2(i) = ms0.prop.rnu2;
    rnu3(i) = ms0.prop.rnu3;
    rnu4(i) = ms0.prop.rnu4;
    aO4C1C2(i) = ms0.prop.aO4C1C2;
    aC1C2C3(i) = ms0.prop.aC1C2C3;
    aC2C3C4(i) = ms0.prop.aC2C3C4;
    aC3C4O4(i) = ms0.prop.aC3C4O4;
    aC4O4C1(i) = ms0.prop.aC4O4C1;
   
    sumofnu(i) = ms0.prop.sumofnu;

    Pdeg(i) = ms0.prop.Pdeg;
    numax(i) = ms0.prop.numax;

    if isfield(ms0.prop,'O5Hn')
        O5Hn(i) = ms0.prop.O5Hn;
    else
        O5Hn(i) = 0;
    end
    if isfield(ms0.prop,'rO5H5')
      rO5H5(i) = ms0.prop.rO5H5;
    else
      rO5H5(i) = 1;
    end
    if isfield(ms0.prop,'rO3H3')
      rO3H3(i) = ms0.prop.rO3H3;
    else
      rO3H3(i) = 0;
    end

    if isfield(ms0.SP,'energy') 
        SPenergy(i) = ms0.SP.energy;
    else
        SPenergy(i) = 0;
    end;
    if isfield(ms0.GO,'energy')
        GOenergy(i) = ms0.GO.energy;
    else
        GOenergy(i) = 0;
    end;
    if ~isempty(ms0.DM), DM(i) = ms0.DM; else DM(i)=0; end;
    
    if strcmp(usedpackage,'Gaussian')
        ZPE(i)=inf;
        TEC(i)=inf;
        EEC(i)=inf;
        GEC(i)=inf;

        if isfield(ms0.gaussian,'T')
          tind = find(ms0.gaussian.T==T);
          if isempty(tind), error(['Desired temperature T=' num2str(T) 'K is not found']), end  
        else  %???
          tind=1;
        end

        if isfield(ms0.gaussian,'ZPE') && ~isempty(ms0.gaussian.ZPE), ZPE(i) = ms0.gaussian.ZPE; end;

        if isfield(ms0.gaussian,'TEC') && ~isempty(ms0.gaussian.TEC), TEC(i) = ms0.gaussian.TEC(tind); end;
        if isfield(ms0.gaussian,'EEC') && ~isempty(ms0.gaussian.EEC), EEC(i) = ms0.gaussian.EEC(tind); end;
        if isfield(ms0.gaussian,'GEC') && ~isempty(ms0.gaussian.GEC), GEC(i) = ms0.gaussian.GEC(tind); end;
    else
        if isfield(ms0,'ZPE') && ~isempty(ms0.ZPE), ZPE(i) = ms0.ZPE; else ZPE(i)=inf; end;
    end
    if isfield(ms0.GO,'steps'), GOsteps(i) = ms0.GO.steps; else GOsteps(i)=0; end;


    if strcmp(usedpackage,'Gaussian')

      %fill gausenergy substructure with SP results with diffrent theories/basises
      for j=1:numel(fdescgausenergy_all)
        fieldnamebuf=fdescgausenergy_all(j);
        if isfield(ms0.gaussian,fieldnamebuf{1});
          gausenergy(j,i) = ms0.gaussian.(fieldnamebuf{1});
        else
          gausenergy(j,i) = 0;
        end
      end

      if isfield(ms0,'freq')
        freq.freq(i,:)=ms0.freq.freq;

        %zhr100927 low freq modes test 
        freq.lowerKT(i)=sum(ms0.freq.freq(ms0.freq.freq<208.53));
        freq.higherKT(i)=sum(ms0.freq.freq(ms0.freq.freq>=208.53));

        Hind(i,:)=strcmpcellar(ms0.labels,'H');
        Hbondpindlabels(i,:)=pind.labels(ms0.pind(Hind(i,:)));
        if (~isfield(ms0.prop,'HbondLMind')) || (max(ms0.prop.HbondLMind)==0) %when no frequency information or it is zeros
          HbondHind(i)={[]};
          HbondLMind(i,:)=zeros(1,size(Hind,2));
          HbondVMind(i,:)=zeros(1,size(Hind,2));
%          HbondLMind(i,:)=1:size(Hind,2);
%          HbondVMind(i,:)=1:size(Hind,2);
        else
          HbondHind(i)={ms0.prop.HbondHind};

          HbondLMind(i,:)=ms0.prop.HbondLMind;
          HbondVMind(i,:)=ms0.prop.HbondVMind;
        end

        if isfield(ms0.freq,'intencity')
          freq.intensity(i,:)=ms0.freq.intencity;
          freq.forceconst(i,:)=ms0.freq.forceconst;
        else
          freq.intensity(i,:)=zeros(1,ms0.atomnum*3-6);
          freq.forceconst(i,:)=zeros(1,ms0.atomnum*3-6);
        end

      end
    end


    %--distance and angles & torsion matrixes
    [distmatrorder2,distmatrordii2] = sortrows(ms0.prop.distmatrpinds);
    if ~all(distmatrorder==distmatrorder2)
    error('Pinds of distance labels are not identical!');
    end 
    distmatr(i,:) = ms0.prop.distmatr(distmatrordii2);

    [anglematrorder2,anglematrordii2] = sortrows(ms0.prop.anglematrpinds);
    if ~all(anglematrorder==anglematrorder2)
    error('Pinds of angle labels are not identical!');
    end 
    anglematr(i,:) = ms0.prop.anglematr(anglematrordii2);

    [tormatrorder2,tormatrordii2] = sortrows(ms0.prop.tormatrpinds);
    if ~all(tormatrorder==tormatrorder2)
    error('Pinds of torsion labels are not identical!');
    end 
    tormatr(i,:) = ms0.prop.tormatr(tormatrordii2);
    %-----------------------------------------------------------------------------


    Hbonddata(i,ismember(Hbonddatadesc,ms0.prop.Hbonddatadesc))=ms0.prop.Hbonddata;
    nbonds(i)=numel(ms0.prop.Hbonddatadesc)/2;  % H bond number for current molecule

end %for i=1:recnum
clear ms0

%if all ZPE or TEc are infinite than change them to 0
if ~sum(isfinite(ZPE))
  ZPE=zeros(size(ZPE));
  warning('rot09:ZPEzeros','ZPE values are zeros!');
end 
if ~sum(isfinite(TEC))
  TEC=zeros(size(TEC));
  warning('rot09:TECzeros','TEC values are zeros!');
end 
if ~sum(isfinite(EEC))
  EEC=zeros(size(EEC));
  warning('rot09:EECzeros','EEC values are zeros!');
end 
if ~sum(isfinite(GEC))
  GEC=zeros(size(GEC));
  warning('rot09:GECzeros','GEC values are zeros!');
end 

if strcmp(usedpackage,'Gaussian')

%[HbondLMfreq,HbondVMfreq]=deal(zeros(size(HbondLMind)));
   for i=1:recnum
     ind=zeros(0);
     if ~((moltype>100 && any(mod(moltype,100)==[30,32])) || moltype==22) %non 5'nucleotides
       iH53= workdb(i).ind(find(find(strcmp(pind.labels,'pH53'))==workdb(i).pind)); 
       ind(end+1)=find(Hind(i,:)==iH53);
     end
     if ~((moltype>100 && any(mod(moltype,100)==[31,12])) || any(moltype==[23,413])) %non 3'nucleotides or dideoxy molecules or AZT
       iH32= workdb(i).ind(find(find(strcmp(pind.labels,'pH32'))==workdb(i).pind)); 
       ind(end+1)=find(Hind(i,:)==iH32);
     end
     if any(moltype==[7,14])
       iH22= workdb(i).ind(find(find(strcmp(pind.labels,'pH22'))==workdb(i).pind)); 
       ind(end+1)=find(Hind(i,:)==iH22);
     end

     
     [Hbondfreqdesc,sortind]=sort(Hbondpindlabels(i,:)); %this is for case when atoms in different structures are in different order
%     Hbondpindlabels(i,:)=Hbondpindlabels(i,sortind);
     Hind(i,:)=Hind(i,sortind);  
     
     if all(HbondLMind(i,:)) %HbondLMind contains no zeros
         HbondLMind(i,:)=HbondLMind(i,sortind);  
         HbondVMind(i,:)=HbondVMind(i,sortind);

         HbondLMfreq(i,:)=freq.freq(i,HbondLMind(i,:));
         HbondVMfreq(i,:)=freq.freq(i,HbondVMind(i,:));
         HbondVMintensity(i,:)=freq.intensity(i,HbondVMind(i,:));

         Hbondindprint(i,:) = [HbondLMind(i,ind) HbondVMind(i,ind)];
         HbondFC(i,:) = freq.forceconst(i,HbondLMind(i,ind));
     else
         Hbondindprint(i,:) = 0;
         HbondFC(i,:) = 0;
     end
   end %i=1:recnum

   Hbondindprintdesc = [{'librH53ind'},{'librH32ind'},{'valH53ind'},{'valH32ind'}];
   HbondFCdesc = [{'librH53 FC,mDyne/A'},{'librH32 FC,mDyne/A'}];
   if any(moltype==[7,14])
     Hbondindprintdesc = [{'librH53ind'},{'librH32ind'},{'librH22ind'},{'valH53ind'},{'valH32ind'},{'valH22ind'}];
     HbondFCdesc = [HbondFCdesc {'librH22 FC,mDyne/A'}];
   end  

   [HbondLMfreqdesc,HbondVMfreqdesc,HbondVMintensitydesc]=deal([]);
   for i=1:numel(Hbondfreqdesc)
     HbondLMfreqdesc=[HbondLMfreqdesc {['libr ' Hbondfreqdesc{i} ' cm^-1']}];
     HbondVMfreqdesc=[HbondVMfreqdesc {['val ' Hbondfreqdesc{i} ' cm^-1']}];
     HbondVMintensitydesc=[HbondVMintensitydesc {['Ival ' Hbondfreqdesc{i} ' KM/mol']}];
   end


   atomsinHbond = zeros(size(Hind));
   for i=1:recnum
     if ~isempty(HbondHind{i})
       [A,B]=meshgrid(Hind(i,:),unique(HbondHind{i})); %unique is for choice when HbondHind{i} contains H atoms involved in several Hbonds
       atomsinHbond(i,:)=sum(A==B,1);
     end
   end

   Hbondmodesnum1=sum(atomsinHbond,1);
   [HbondLMfreq1,HbondVMfreq1,HbondVMintensity1]=deal(zeros(size(Hind)));
   inds=find(atomsinHbond);
   HbondLMfreq1(inds)=HbondLMfreq(inds);
   HbondVMfreq1(inds)=HbondVMfreq(inds);
   HbondVMintensity1(inds)=HbondVMintensity(inds);

   warning off MATLAB:divideByZero
   meanLM1 = sum(HbondLMfreq1,1)./Hbondmodesnum1; %mean value of libration modes frequinces included in H bonds
   meanVM1 = sum(HbondVMfreq1,1)./Hbondmodesnum1; %mean value of stretching(valency) modes frequinces included in H bonds
   
   Hbondmodesnum0=recnum-Hbondmodesnum1;
   HbondLMfreq0=HbondLMfreq-HbondLMfreq1;
   HbondVMfreq0=HbondVMfreq-HbondVMfreq1;
   HbondVMintensity0=HbondVMintensity-HbondVMintensity1;

   meanLM0 = sum(HbondLMfreq0,1)./Hbondmodesnum0; %mean value of libration modes frequinces not included in H bonds
%   meanVM0 = sum(HbondVMfreq0,1)./Hbondmodesnum0; %mean value of stretching(valency) modes frequinces not included in H bonds
   meanVM0 = max(HbondVMfreq0,[],1);
   HbondVMintensity0(HbondVMintensity0==0)=NaN; %do not take into account zero values of intensity
   meanVMint0 = min(HbondVMintensity0,[],1);
   warning on MATLAB:divideByZero
   
   dmeanLM = meanLM1-meanLM0;
   dmeanVM = meanVM1-meanVM0;


   gausenergy(gausenergy==0)=inf;

   iconfdesc=1;
   ibonddesc=2;
   irOO=3;
   irHO=4;
   irOH=5;
   idrOH=6;
   iaOHO=7;
   ifreqVM=8;
   idfreqVM=9;
   ifreqLM=10;
   iintVM=11;
   irelintVM=12; % Ivm/Ivm0

   tableHbonddesc=[{'confdesc'},{'Hbonddesc'},{'rAB'},{'rHB'},{'rAH'},{'rdAH'},{'aAHB'},{'freqVM'},{'dfreqVM'},...
                   {'freqLM'},{'intVM'},{'intVM/intVM0'}];

   bondpinds=zeros(0,2);

   tableHbond={};
   for i=1:recnum
    N=numel(workdb(i).prop.HbondHind);
      for j=1:N

        bondHpind=workdb(i).pind(workdb(i).prop.HbondO1ind(j));
        bondOpind=workdb(i).pind(workdb(i).prop.HbondHind(j));
        bondpinds(end+1,:)=[bondOpind bondHpind];

        if all(workdb(i).labels{workdb(i).prop.HbondO2ind(j)}~='C') && all(workdb(i).labels{workdb(i).prop.HbondO1ind(j)}~='C')
            tp=1; %table part: 1 - is for ordinary H-bonds, 2 - for non
        else
            tp=2;
        end 
        
        itableHbond(tp)=itableHbond(tp)+1;
        tableHbond(iconfdesc,itableHbond(tp),tp)={workdb(i).prop.sdesc};
        tableHbond(ibonddesc,itableHbond(tp),tp)=workdb(i).prop.Hbonddatadesc(j);
        tableHbond(irOO,itableHbond(tp),tp)={workdb(i).prop.Hbonddata(N+j)};
    
        indH=find(Hind(i,:)==workdb(i).prop.HbondHind(j)); %position of current H atom's data value in Hbond* vectors

        tableHbond(irHO,itableHbond(tp),tp)={adist(workdb(i),workdb(i).prop.HbondHind(j),workdb(i).prop.HbondO2ind(j))};
        tableHbond(irOH,itableHbond(tp),tp)={adist(workdb(i),workdb(i).prop.HbondO1ind(j),workdb(i).prop.HbondHind(j))};
        tableHbond(iaOHO,itableHbond(tp),tp)={workdb(i).prop.Hbonddata(j)};
    
        if HbondVMind(i,indH)
            tableHbond(ifreqVM,itableHbond(tp),tp)={freq.freq(i,HbondVMind(i,indH))};
            tableHbond(ifreqLM,itableHbond(tp),tp)={freq.freq(i,HbondLMind(i,indH))};
            tableHbond(iintVM,itableHbond(tp),tp)={freq.intensity(i,HbondVMind(i,indH))};
        else
            tableHbond(ifreqVM,itableHbond(tp),tp)={0};
            tableHbond(ifreqLM,itableHbond(tp),tp)={0};
            tableHbond(iintVM,itableHbond(tp),tp)={0};
        end
        tableHbond(irelintVM,itableHbond(tp),tp)={0};


     end
   end  
   OHinbondsind=unique(bondpinds,'rows'); %OH atoms pind indexes that are involved in Hbonds
   OHinbondsnum=size(OHinbondsind,1);
   rOHall=[];
   for i=1:recnum
   for j=1:OHinbondsnum
     rOHall(i,j)=adist(workdb(i),find(workdb(i).pind==OHinbondsind(j,1)),find(workdb(i).pind==OHinbondsind(j,2)));
%    freqVMall(i,j)=ms0.
   end
   end
   OHmindist=min(rOHall);
   OHdistdelta=max(rOHall)-min(rOHall);
   [rOHallsorted,rOHallsortedind]=sort(rOHall);

   OHminmeandist=0;
   if ~isempty(rOHall)
       elnumformeansearch=uint16(min(max(0.15*recnum,10),recnum));  %select number of least OH distance values to calculate mean value
       OHminmeandist=mean(rOHallsorted(1:elnumformeansearch,:));
   end

   if recnum>1
       itableHbond=[0 0];
    %   tableHbond=[];
       for i=1:recnum
          N=numel(workdb(i).prop.HbondHind);
          for j=1:N
            if all(workdb(i).labels{workdb(i).prop.HbondO2ind(j)}~='C') && all(workdb(i).labels{workdb(i).prop.HbondO1ind(j)}~='C')
                tp=1; %table part: 1 - is for ordinary H-bonds, 2 - for non
            else
                tp=2;
            end 
            itableHbond(tp)=itableHbond(tp)+1;
            OHbond=repmat([workdb(i).pind(workdb(i).prop.HbondHind(j)) workdb(i).pind(workdb(i).prop.HbondO1ind(j))],OHinbondsnum,1);
            OHbondind=find(sum(OHinbondsind==OHbond,2)==2);
            tableHbond(idrOH,itableHbond(tp),tp)={tableHbond{irOH,itableHbond(tp),tp}-OHminmeandist(OHbondind)};

         end
       end  
   end
       
%   if N
   if itableHbond(2)
       tableHbond=[tableHbond(:,1:itableHbond(1),1) tableHbond(:,1:itableHbond(2),2)]; %rearrange table: at first ordinary Hbonds then others
   end
%   end

end


fdesc00=[{'#'},{'desc'},{'sdesc'},{'sdescnew'},{'fulldesc'},{'new'},...
     {'filename'},{'worktitle'},{'date'},{'comment'}];
fdesc0=[{'hta2'},{'hta3'},{'hta4'},{'hta5'},{'hta0'}];
fdesc=[{'tbeta'},{'tgamma'},{'tdelta'},{'tepsilon'},{'teta'},...
    {'tteta'},{'tchi'},{'tkappa2'},...
    {'tnu0'},{'tnu1'},{'tnu2'},{'tnu3'},{'tnu4'},{'sumofnu'},{'Pdeg'},{'numax'},...
    {'rnu0,A'},{'rnu1,A'},{'rnu2,A'},{'rnu3,A'},{'rnu4,A'},...
    {'drnu0,%'},{'drnu1,%'},{'drnu2,%'},{'drnu3,%'},{'drnu4,%'},...
    {'rO3H3, A'},{'rO5H5,A'},{'C2 dist to ring plane, A'},{'C3 dist to ring plane, A'},{'rglyc'},...
    {'tau10,deg'},{'tau01,deg'},{'tau02,deg'},{'tau03,deg'},{'max tor ang in base,deg'},...
    torbpldesc,...
    {'tau11,C5C4N4H41'},{'tau12,C5C4N4H42'}];
%    {'C6N1C2N3'},{'N1C2N3C4'},{'C2N3C4C5'},{'N3C4C5C6'},{'C4C5C6N1'},{'C5C6N1C2'},...
fdesc_=[{'H5 C4'},{'O5 C3'},{'O4 O3'},{'C4 H3'},{'C3 H2'}...
    {'O3 O2'},{'O4''C2'},{'H31H32'},...
    {'C4 C2'},{'O4 C3'},{'C1 C4'},{'C2 O4'},{'C3 C1'},{''},{''},{''},...
    {'O4 C1'},{'C1 C2'},{'C2 C3'},{'C3 C4'},{'C4 O4'}];
fdescen=[{'GO energy, a.e.'},{'SP energy, a.e.'},{'Dipole moment, Debay'},...
     {'ZPE, kcal/mol'},{'TEC, kcal/mol'},{'EEC, kcal/mol'},{'GEC, kcal/mol'},{'GO steps #'}];



indmain=['A3';'B3';'C3';'D3';'E3';'F3';'G3';'H3';'I3';'J3';'K3';...
         'L3';'M3';'N3';'O3';'P3';'Q3';'R3';'S3';'T3';'U3';'V3';...
         'W3';'X3';'Y3';'Z3'];
indgeom=['AA3';'AB3';'AC3';'AD3';'AE3';'AF3';'AG3';'AH3';'AI3';'AJ3';...
         'AK3';'AL3';'AM3';'AN3';'AO3';'AP3';'AQ3';'AR3';'AS3';'AT3';...
         'AU3';'AV3';'AW3';'AX3';'AY3';'AZ3';'BA3';'BB3';'BC3';'BD3';...
         'BE3';'BF3';'BG3';'BH3';'BI3';'BJ3';...
         'BK3';'BL3';'BM3';'BN3';'BO3';'BP3';'BQ3';'BR3';'BS3';'BT3';...
         'BU3';'BV3';'BW3';'BX3';'BY3';'BZ3'];
inden  =['CA3';'CB3';'CC3';'CD3';'CE3';'CF3';'CG3';'CH3';'CI3';'CJ3';...
         'CK3';'CL3';'CM3';'CN3';'CO3';'CP3';'CQ3';'CR3';'CS3';'CT3';...
         'CU3';'CV3';'CW3';'CX3';'CY3';'CZ3'];
indcalc=['DA3';'DB3';'DC3';'DD3';'DE3';'DF3';'DG3';'DH3';'DI3';'DJ3';...
         'DK3';'DL3';'DM3';'DN3';'DO3';'DP3';'DQ3';'DR3';'DS3';'DT3';...
         'DU3';'DV3';'DW3';'DX3';'DY3';'DZ3'];
indE   =['EA3';'EB3';'EC3';'ED3';'EE3';'EF3';'EG3';'EH3';'EI3';'EJ3';...
         'EK3';'EL3';'EM3';'EN3';'EO3';'EP3';'EQ3';'ER3';'ES3';'ET3';...
         'EU3';'EV3';'EW3';'EX3';'EY3';'EZ3'];
indHbond=['FA3';'FB3';'FC3';'FD3';'FE3';'FF3';'FG3';'FH3';'FI3';'FJ3';...
          'FK3';'FL3';'FM3';'FN3';'FO3';'FP3';'FQ3';'FR3';'FS3';'FT3';...
          'FU3';'FV3';'FW3';'FX3';'FY3';'FZ3';...
          'GA3';'GB3';'GC3';'GD3';'GE3';'GF3';'GG3';'GH3';'GI3';'GJ3';...
          'GK3';'GL3';'GM3';'GN3';'GO3';'GP3';'GQ3';'GR3';'GS3';'GT3';...
          'GU3';'GV3';'GW3';'GX3';'GY3';'GZ3'];


if flwritefile %write file only if flag set
  if exist(xlsfile,'file')==2
      dlm=strfind(xlsfile,'.');
      xlsfileold=[xlsfile(1:dlm(end)-1) '~' xlsfile(dlm(end):end)];
      copyfile(xlsfile,xlsfileold);
  end

  Excel = actxserver ('Excel.Application');
  File=fullfile(pwd,xlsfile);
  if ~exist(File,'file')
    ExcelWorkbook = Excel.workbooks.Add;
    ExcelWorkbook.SaveAs(File,1);
    ExcelWorkbook.Close(false);
  end
  invoke(Excel.Workbooks,'Open',File);


  xlswrite1spec(xlsfile,fdesc00,'ribose','A1');
  xlswrite1spec(xlsfile,[fdesc0 fdesc],'ribose','AA1');
  xlswrite1spec(xlsfile,[{''} {''} {''} {''} {'C5C2'} fdesc_],'ribose','AA2');
  xlswrite1spec(xlsfile,fdescen,'ribose','CA1');

  i=1;
  xlswrite1spec(xlsfile,(1:recnum)','ribose',indmain(i,:)); i=i+1;
  descdata=[desc; sdesc; sdescnew; confdesc; new; filename; worktitle; ldate];
  xlswrite1spec(xlsfile,descdata','ribose',indmain(i,:)); i=i+size(descdata,1);

  maxcommentsize=800;
  for comind=1:numel(comment)
      comsize=numel(comment{comind});
      if comsize>maxcommentsize
          comment(comind)={comment{comind}(1:maxcommentsize)};
          warning(['Comment of record #' int2str(comind) ' is cut off from ' int2str(comsize) ' to ' int2str(maxcommentsize)]);
      end
  end
  xlswrite1spec(xlsfile,comment','ribose',indmain(i,:));

  rnu=[rnu0;rnu1;rnu2;rnu3;rnu4];
  mean_rnu=mean(rnu,2);
  min_rnu=min(rnu,[],2);
  max_rnu=max(rnu,[],2);
  d_rnu=(rnu-repmat(mean_rnu,1,recnum))./rnu*100; %relative
  std_rnu=std(rnu,1,2); %sqrt(mean((x-mean(x)).^2))
  relstd_rnu=std_rnu./mean_rnu;

  aalfa=[aO4C1C2;aC1C2C3;aC2C3C4;aC3C4O4;aC4O4C1];
  mean_aalfa=mean(aalfa,2);
  min_aalfa=min(aalfa,[],2);
  max_aalfa=max(aalfa,[],2);
  d_aalfa=(aalfa-repmat(mean_aalfa,1,recnum))./aalfa*100;
  std_aalfa=std(aalfa,1,2); %sqrt(mean((x-mean(x)).^2))
  relstd_aalfa=std_aalfa./mean_aalfa;

  tcommonvalues={};
  tcommonvalues(1,2:6)=[{'rnu0'},{'rnu1'},{'rnu2'},{'rnu3'},{'rnu4'}];
  tcommonvalues(2:6,1)=[{'min,A'},{'max,A'},{'mean,A'},{'std,A'},{'relstd'}];
  tcommonvalues(2:6,2:6)=num2cell([min_rnu,max_rnu,mean_rnu,std_rnu,relstd_rnu]');
  tcommonvalues(7,:)=cell(1,6);
  tcommonvalues(8,2:6)=[{'aO4C1C2'},{'aC1C2C3'},{'aC2C3C4'},{'aC3C4O4'},{'aC4O4O1'}];
  tcommonvalues(9:13,1)=[{'min,deg'},{'max,deg'},{'mean,deg'},{'std,deg'},{'relstd'}];
  tcommonvalues(9:13,2:6)=num2cell([min_aalfa,max_aalfa,mean_aalfa,std_aalfa,relstd_aalfa]');

  if strcmp(usedpackage,'Gaussian')

    [R,C]=size(tcommonvalues);
    tcommonvalues(R+1,:)=cell(1,C);
    [R,C]=size(tcommonvalues);
    desc=[{'mean libration modes freq (not in H bonds)'},HbondLMfreqdesc,...
          {'mean libration modes freq (in H bonds)'},HbondLMfreqdesc,...
          {'difference between mean libration frequences'},HbondLMfreqdesc,...
          {'mean stretching modes freq (not in H bonds)'},HbondVMfreqdesc,...
          {'mean stretching modes freq (in H bonds)'},HbondVMfreqdesc,...
          {'difference between mean stretching frequences'},HbondVMfreqdesc,...
          {'minimal intensity of modes not in H bonds, KM/Mole'},HbondVMintensitydesc]';
    data=[NaN,meanLM0,NaN,meanLM1,NaN,dmeanLM,NaN,meanVM0,NaN,meanVM1,NaN,dmeanVM,NaN,meanVMint0]';
    tcommonvalues(R+1:R+numel(desc),1:2)=[desc num2cell(data)];

    xlswrite1spec(xlsfile,[tableHbonddesc;tableHbond'],'Hbondtable','A1');

  end
  
  if 1  % write statistics over torsion sectors
      [R,C]=size(tcommonvalues);
      tcommonvalues(R+1,:)=cell(1,C);
      R=R+1;
      tcommonvalues(R+1,2:5)=[{'min'},{'max'},{'mean'},{'std'}];
      tcommonvalues((R+2):(R+17),1)=[{'syn'},{'anti'},{'north'},{'south'},...
                                    {'betagp'},{'betat'},{'betagm'},...
                                    {'gammagp'},{'gammat'},{'gammagm'},...
                                    {'deltagp'},{'deltat'},{'deltagm'},...
                                    {'epsilongp'},{'epsilont'},{'epsilongm'},...
                                    ];

      val_syn=[];
      val_anti=[];
      if ~any(moltype==[7,8])
          for ind=1:numel(sdesc)
              ttt=sdesc{ind}(end);
              if ttt=='S' || ttt=='T' || ttt=='V'
                 val_syn(end+1)=tchi(ind);
              elseif ttt=='A' || ttt=='B' || ttt=='C'
                 val_anti(end+1)=tchi(ind);
              else
                 error(['Unknown conformation chintype detected: ' workdb(curconf).prop.sdesc]);
              end
          end
      end
      val=val_syn;
      if isempty(val), val=0;, end;
      syn_vals = [min(val) max(val) mean(val) std(val)];
      val=val_anti;
      if isempty(val), val=0;, end;
      anti_vals = [min(val) max(val) mean(val) std(val)];

      val = Pdeg(find(Pdeg<90 | Pdeg>270)); %P in [0,360]
      val(val>180)=val(val>180)-360;
      if isempty(val), val=0;, end;
      north_vals = [min(val) max(val) mean(val) std(val)];
      val = Pdeg(find(Pdeg>=90 & Pdeg<=270));
      if isempty(val), val=0;, end;
      south_vals = [min(val) max(val) mean(val) std(val)];

      %beta in [-180,180]
      val = tbeta(find(ceil(tbeta/120)==1)); % (0,120]
      if isempty(val), val=0;, end;
      tbetagp_vals = [min(val) max(val) mean(val) std(val)];
      val = tbeta(find(ceil(tbeta/120)==-1 ))+360; % (-240,-120]
      val = [val tbeta(find(ceil(tbeta/120)==2))]; % (120,240] 
      if isempty(val), val=0;, end;
      tbetat_vals = [min(val) max(val) mean(val) std(val)];
      val = tbeta(find(ceil(tbeta/120)==0)); % (-120,0]
      if isempty(val), val=0;, end;
      tbetagm_vals = [min(val) max(val) mean(val) std(val)];

      val = tgamma(find(ceil(tgamma/120)==1));
      if isempty(val), val=0;, end;
      tgammagp_vals = [min(val) max(val) mean(val) std(val)];
%      val = tgamma(find(ceil(tgamma/120)==2 | ceil(tgamma/120)==-1));
      val = tgamma(find(ceil(tgamma/120)==-1 ))+360; % (-240,-120]
      val = [val tgamma(find(ceil(tgamma/120)==2))]; % (120,240] 
      if isempty(val), val=0;, end;
      tgammat_vals = [min(val) max(val) mean(val) std(val)];
      val = tgamma(find(ceil(tgamma/120)==0));
      if isempty(val), val=0;, end;
      tgammagm_vals = [min(val) max(val) mean(val) std(val)];

      val = tdelta(find(ceil(tdelta/120)==1));
      if isempty(val), val=0;, end;
      tdeltagp_vals = [min(val) max(val) mean(val) std(val)];
%      val = tdelta(find(ceil(tdelta/120)==2 | ceil(tdelta/120)==-1));
      val = tdelta(find(ceil(tdelta/120)==-1 ))+360; % (-240,-120]
      val = [val tdelta(find(ceil(tdelta/120)==2))]; % (120,240] 
      if isempty(val), val=0;, end;
      tdeltat_vals = [min(val) max(val) mean(val) std(val)];
      val = tdelta(find(ceil(tdelta/120)==0));
      if isempty(val), val=0;, end;
      tdeltagm_vals = [min(val) max(val) mean(val) std(val)];

      val = tepsilon(find(ceil(tepsilon/120)==1));
      if isempty(val), val=0;, end;
      tepsilongp_vals = [min(val) max(val) mean(val) std(val)];
%      val = tepsilon(find(ceil(tepsilon/120)==2 | ceil(tepsilon/120)==-1));
      val = tepsilon(find(ceil(tepsilon/120)==-1 ))+360; % (-240,-120]
      val = [val tepsilon(find(ceil(tepsilon/120)==2))]; % (120,240] 
      if isempty(val), val=0;, end;
      tepsilont_vals = [min(val) max(val) mean(val) std(val)];
      val = tepsilon(find(ceil(tepsilon/120)==0));
      if isempty(val), val=0;, end;
      tepsilongm_vals = [min(val) max(val) mean(val) std(val)];

      tcommonvalues((R+2):(R+17),2:5)=num2cell([
        syn_vals; anti_vals; north_vals; south_vals;...
        tbetagp_vals; tbetat_vals; tbetagm_vals;...
        tgammagp_vals; tgammat_vals; tgammagm_vals;...
        tdeltagp_vals; tdeltat_vals; tdeltagm_vals;...
        tepsilongp_vals; tepsilont_vals; tepsilongm_vals;...
        ]);
  end
  xlswrite1spec(xlsfile,tcommonvalues,'commonvalues','A1');

  tdistmatr=[{'sdesc'} distmatrdesc'; sdesc' num2cell(distmatr)];
  tanglematr=[anglematrdesc' ; num2cell(anglematr) ];
  ttormatr=[tormatrdesc' ; num2cell(tormatr) ];
  xlswrite1spec(xlsfile,[tdistmatr tanglematr ttormatr [{'sdesc'}; sdesc']],'dist&angle matrix','A1');

  %---distance & angles & torsions correlations
  corrdata=[distmatr anglematr tormatr];
  corrdatadesc=[distmatrdesc; anglematrdesc; tormatrdesc];
  if size(corrdata,1)==1
      [r,p] = deal(zeros(size(corrdata,2)));
  else
      [r,p] = corrcoef(corrdata);
  end
  %[i,j] = find(p<0.00005); %matrix of p-values for testing
  %    the hypothesis of no correlation.  Each p-value is the probability
  %    of getting a correlation as large as the observed value by random
  %    chance, when the true correlation is zero.
%  r(p>0.0000001)=NaN;

%  r(p>0.1)=NaN;
  ii=find(abs(r)<0.05);
  r(r==1)=NaN;
  rr=num2cell(r);
  rr(ii)={'-'};
  rr=reshape(rr,size(r));

  ttt=repmat({''},1,size(rr,1));
  tcorrcoefmatr=[corrdatadesc ttt' rr];
  tcorrcoefmatr=[{''} {''} corrdatadesc'; {''} {''} ttt; tcorrcoefmatr];
  xlswrite1spec(xlsfile,tcorrcoefmatr,'geomcorr matrix','A1');

  %---torsions curcular correlations 
  rowsused=size(rr,1)+2+1+1;

  corrdata=tormatr/180*pi; %to radians
  corrdatadesc=tormatrdesc;
  r = corrcoefcirc(corrdata);
  r(abs(r)<0.7)=NaN;
  r(r==1)=NaN;
  rr=num2cell(r);
  rr=reshape(rr,size(r));

  ttt=repmat({''},1,size(rr,1));
  tcorrcoefmatr=[corrdatadesc ttt' rr];
  tcorrcoefmatr=[{''} {''} corrdatadesc'; {''} {''} ttt; tcorrcoefmatr];
  xlswrite1spec(xlsfile,tcorrcoefmatr,'geomcorr matrix',['A' int2str(rowsused)]);

  %---main torsions correlation begin
  rowsused=rowsused+size(tcorrcoefmatr,1)+2;

  corrdatadesc=[{'nu0'}; {'nu1'}; {'nu2'}; {'nu3'}; {'nu4'}; {'gamma'}; {'epsilon'}; {'beta'}; {'delta'}];
  if ~any(moltype==[7,8])
    corrdatadesc=[corrdatadesc; {'chi'}];
  end
  torsigns = ones(size(corrdatadesc));
  maintorinds = zeros(size(corrdatadesc));
  
  if mod(moltype,100)==11 %tionucleosides
      maintorinds(1)=strcmpcellar(tormatrdesc,'pC2pC1pS4pC4'); %nu0
      maintorinds(2)=strcmpcellar(tormatrdesc,'pC3pC2pC1pS4'); %nu1
      maintorinds(3)=strcmpcellar(tormatrdesc,'pC1pC2pC3pC4'); %nu2
      maintorinds(4)=strcmpcellar(tormatrdesc,'pC2pC3pC4pS4'); %nu3
      maintorinds(5)=strcmpcellar(tormatrdesc,'pC1pS4pC4pC3'); %nu4
  else
      maintorinds(1)=strcmpcellar(tormatrdesc,'pC2pC1pO4pC4'); %nu0
      maintorinds(2)=strcmpcellar(tormatrdesc,'pC3pC2pC1pO4'); %nu1
      maintorinds(3)=strcmpcellar(tormatrdesc,'pC1pC2pC3pC4'); %nu2
      maintorinds(4)=strcmpcellar(tormatrdesc,'pC2pC3pC4pO4'); %nu3
      maintorinds(5)=strcmpcellar(tormatrdesc,'pC1pO4pC4pC3'); %nu4
  end

  torsigns(1)=-1;
  torsigns(2)=-1;
  torsigns(5)=-1;
  maintorinds(6)=strcmpcellar(tormatrdesc,'pC3pC4pC5pO5'); %gamma
  torsigns(6)=-1;
  
  if mod(moltype,100)==12 %dideoxy molecules
      maintorinds(7)=1; %need to be corrected !!
  elseif ((moltype>100 && mod(moltype,100)==31) || moltype==23) %3'nucleotides
      maintorinds(7)=strcmpcellar(tormatrdesc,'pC4pC3pO3pP3'); %epsilon
  elseif (moltype==413) %AZT
      maintorinds(7)=strcmpcellar(tormatrdesc,'pC4pC3pN31pN32'); %epsilon
  else
      maintorinds(7)=strcmpcellar(tormatrdesc,'pC4pC3pO3pH32'); %epsilon
  end

  if ((moltype>100 && any(mod(moltype,100)==[30,32])) || moltype==22) %5'nucleotides
      maintorinds(8)=strcmpcellar(tormatrdesc,'pC4pC5pO5pP5'); %beta
  else
      maintorinds(8)=strcmpcellar(tormatrdesc,'pC4pC5pO5pH53'); %beta
  end
  torsigns(8)=-1;

  if mod(moltype,100)==12 %dideoxy molecules
      maintorinds(9)=1; %need to be corrected !!
  elseif mod(moltype,100)==11 %tionucleosides
      maintorinds(9)=strcmpcellar(tormatrdesc,'pO3pC3pC4pS4'); %delta
      torsigns(9)=-1;
  elseif moltype==413 %AZT
      maintorinds(9)=strcmpcellar(tormatrdesc,'pO4pC4pC3pN31'); %delta
  else
      maintorinds(9)=strcmpcellar(tormatrdesc,'pO4pC4pC3pO3'); %delta
  end
    
  if any(moltype==[11,9,12,13,14,21,220,420,520,17,18,19,21]) ||...
     (moltype>=200 && moltype<=299) || (moltype>=400 && moltype<=499) || (moltype>=500 && moltype<=599)
    if mod(moltype,100)==11 %tionucleosides
      maintorinds(10)=strcmpcellar(tormatrdesc,'bC2bN1pC1pS4'); %chi
      torsigns(10)=-1;
    else
      maintorinds(10)=strcmpcellar(tormatrdesc,'pO4pC1bN1bC2'); %chi
    end
  elseif any(moltype==[15,16,110,310]) || (moltype>=100 && moltype<=199) || (moltype>=300 && moltype<=399)
    if mod(moltype,100)==11 %tionucleosides
      maintorinds(10)=strcmpcellar(tormatrdesc,'bC4bN9pC1pS4'); %chi
      torsigns(10)=-1;
    else
      maintorinds(10)=strcmpcellar(tormatrdesc,'pO4pC1bN9bC4'); %chi
    end
  elseif any(moltype==[7,8])
    %do nothing
  else
    warning('rot09:unknownstructure','Unknown structure: Dont know how to calculate chi torsion :(');
  end
  if any(maintorinds==0)
      error('Main torsions correlation: zero indexes detected');
  end
  
  corrdata=tormatr(:,maintorinds)/180*pi; %to radians
  torsignsrep=repmat(torsigns',size(corrdata,1),1);
  corrdata=corrdata.*torsignsrep;
  r = corrcoefcirc(corrdata);
%  r(find(abs(r)<0.5))=NaN;
%  ii=find(abs(r)<0.7);
  rr=num2cell(r);
  %rr(ii)={'-'};
  rr=reshape(rr,size(r));

  ttt=repmat({''},1,size(rr,1));
  tcorrcoefmatr=[corrdatadesc ttt' rr];
  tcorrcoefmatr=[{''} {''} corrdatadesc'; {''} {''} ttt; tcorrcoefmatr];
  xlswrite1spec(xlsfile,tcorrcoefmatr,'geomcorr matrix',['A' int2str(rowsused)]);
  %---main torsions correlation end


%---main torsions correlation - do only for anti conformers begin
  rowsused=rowsused+size(tcorrcoefmatr,1)+2;

  corrdata=tormatr(abs(tchi)>=90,maintorinds)/180*pi; %to radians
  if ~isempty(corrdata)
      torsignsrep=repmat(torsigns',size(corrdata,1),1);
      corrdata=corrdata.*torsignsrep;
      r = corrcoefcirc(corrdata);
  else
      warning('rot09:noanti','There are no ANTI conformers');
      r=zeros(size(corrdata,2));
  end
  rr=num2cell(r);
  rr=reshape(rr,size(r));
  ttt=repmat({''},1,size(rr,1));
  tcorrcoefmatr=[corrdatadesc ttt' rr];
  tcorrcoefmatrA=[{''} {''} corrdatadesc'; {''} {''} ttt; tcorrcoefmatr];
%---main torsions correlation - do only for anti conformers end

%---main torsions correlation - do only for syn conformers begin
  corrdata=tormatr(abs(tchi)<90,maintorinds)/180*pi; %to radians
  if ~isempty(corrdata)
      torsignsrep=repmat(torsigns',size(corrdata,1),1);
      corrdata=corrdata.*torsignsrep;
      r = corrcoefcirc(corrdata);
  else
      warning('rot09:nosyn','There are no SYN conformers');
      r=zeros(size(corrdata,2));
  end
  rr=num2cell(r);
  rr=reshape(rr,size(r));
  ttt=repmat({''},1,size(rr,1));
  tcorrcoefmatr=[corrdatadesc ttt' rr];
  tcorrcoefmatrS=[{''} {''} corrdatadesc'; {''} {''} ttt; tcorrcoefmatr];
   
  xlswrite1spec(xlsfile,[{'for anti conformers:'} ttt {''} {'for syn conformers:'}],'geomcorr matrix',['A' int2str(rowsused-1)]);
  xlswrite1spec(xlsfile,[tcorrcoefmatrA tcorrcoefmatrS],'geomcorr matrix',['A' int2str(rowsused)]);
%---main torsions correlation - do only for syn conformers begin

  
  %---main aminofragments correlations for rCyd begin
  rowsused=rowsused+size(tcorrcoefmatr,1)+3;

  if moltype==240

  corrdatafull=[distmatr anglematr ];
  corrdatadescfull=[distmatrdesc; anglematrdesc];

    corrdatadesc=[{'C4N4'}; {'N4H41'}; {'N4H42'}; {'H41N4H42'}];
    maintorinds = zeros(size(corrdatadesc));
    maintorinds(1)=strcmpcellar(corrdatadescfull,'bC4bN'); %CN
    maintorinds(2)=strcmpcellar(corrdatadescfull,'bNbH41'); %NH1
    maintorinds(3)=strcmpcellar(corrdatadescfull,'bNbH42'); %NH2
    maintorinds(4)=strcmpcellar(corrdatadescfull,'bH41bNbH42'); %H1NH2
    if any(maintorinds==0)
        error('Main torsions correlation: zero indexes detected');
    end
    corrdata=corrdatafull(:,maintorinds)/180*pi; %to radians

    corrdatadesc=[corrdatadesc; {'tau02'}; {'tau03'}];
    corrdata=[corrdata tau02' tau03'];

    r = corrcoefcirc(corrdata);
  %  r(find(abs(r)<0.5))=NaN;
  %  ii=find(abs(r)<0.7);
    rr=num2cell(r);
    %rr(ii)={'-'};
    rr=reshape(rr,size(r));

    ttt=repmat({''},1,size(rr,1));
    tcorrcoefmatr=[corrdatadesc ttt' rr];
    tcorrcoefmatr=[{''} {''} corrdatadesc'; {''} {''} ttt; tcorrcoefmatr];
    xlswrite1spec(xlsfile,tcorrcoefmatr,'geomcorr matrix',['A' int2str(rowsused)]);
  end
  %---main aminofragments correlations for dAdo end

      datageom=[hta2;hta3;hta4;hta5; hta0;...
        tbeta;tgamma;tdelta;tepsilon;teta;tteta;tchi;tkappa2;...
        tnu0;tnu1;tnu2;tnu3;tnu4;sumofnu;Pdeg;numax;rnu0;rnu1;rnu2;rnu3;rnu4;d_rnu;rO3H3;rO5H5;dplC2;dplC3;...
        rglyc;tau10;tau01;tau02;tau03;maxtorinbase;torbpl';tau11;tau12];
      xlswrite1spec(xlsfile,datageom','ribose','AA3');
      
  if all(SPenergy==0) && exist('mostaccenergyfnind','var') && mostaccenergyfnind
    SPenergy=gausenergy(mostaccenergyfnind,:);
  end
  dataen=[GOenergy; SPenergy; DM; ZPE; TEC; EEC; GEC; GOsteps]; %ZPE & TEC are not scaled
  xlswrite1spec(xlsfile,dataen','ribose',inden(1,:));

  %iexclude = find(hta0<150);  %exclude beta-ribose conformations (150 degree is critical value)
  iexclude = strcmpcellar(new,'B');
  GOenergy_ = remels(GOenergy,iexclude);
  SPenergy_ = remels(SPenergy,iexclude);

  ribglvalues.minGOenergy=min(GOenergy_);
  ribglvalues.minSPenergy=min(SPenergy_);
  %ribglvalues.minZPE=min(ZPE);

  fdesccalc=[{'dE GO 6-31G**, kcal/mol'},{'dE SP 6-311++G**, kcal/mol'},...
       {'d(E GO + GEC), kcal/mol'},{'d(E SP + GEC), kcal/mol'}];
  xlswrite1spec(xlsfile,fdesccalc,'ribose','DA1');

  dGOenergykm = CC.encoef*(GOenergy-ribglvalues.minGOenergy); %kcal/mol
  dSPenergykm = CC.encoef*(SPenergy-ribglvalues.minSPenergy); %kcal/mol
  GOenergyGECkm = dGOenergykm+GEC;
  SPenergyGECkm = dSPenergykm+GEC;
  dGOenergyGECkm = GOenergyGECkm-min(GOenergyGECkm);
  dSPenergyGECkm = SPenergyGECkm-min(SPenergyGECkm);


  i=1;
  datacalc=[dGOenergykm; dSPenergykm; dGOenergyGECkm; dSPenergyGECkm];
  xlswrite1spec(xlsfile,datacalc','ribose',indcalc(i,:));i=i+4;
  i=i+1;


  if exist('mostaccenergyfnind','var')
      enfields=[fdescgausenergy_all{mostaccenergyfnind} ', Ee'];
      fdescgausenergy=[fdescgausenergyGECcorr,{enfields},...
                   {'lowest freq mode'},{'lowest mode FC'}];
      if strcmp(usedpackage,'Gaussian')
          xlswrite1spec(xlsfile,fdescgausenergy,'ribose','EA1');

      %    gausenergy_ = remels(gausenergy,iexclude); %exclude beta-ribose conformations (150 degree is critical value)
          mingausenergy=repmat(min(gausenergy,[],2),1,recnum);
          gausenergykm = CC.encoef*(gausenergy-mingausenergy); %kcal/mol
          gausenergyGECkm = gausenergykm+repmat(GEC,size(gausenergy,1),1);  %GEC correcting
          dgausenergyGECkm = gausenergyGECkm-repmat(min(gausenergyGECkm,[],2),1,recnum); %GEC correcting

          xlswrite1spec(xlsfile,[dgausenergyGECkm; gausenergykm(mostaccenergyfnind,:)]','ribose','EA3');
      end
  end

  if strcmp(usedpackage,'Gaussian')
      if isfield(freq,'freq')
        xlswrite1spec(xlsfile,[freq.freq(:,1) freq.forceconst(:,1)],'ribose',indE(numel(fdescgausenergyGECcorr)+2,:)); %lowest frequency mode

        %zhr100927 low freq modes test 
        %sum of modes freqs with energy less then kT=0.6kcal/mol and with energies higher than kT
%        xlswrite1spec(xlsfile,[freq.lowerKT freq.higherKT]','ribose',indE(numel(fdescgausenergyGECcorr)+4,:)); 

      end
  end

  Hbonddata(Hbonddata==0)=NaN; %empties data when there is no H bond
  Hbond_desc1=[{'nbonds'} Hbonddatadesc HbondLMfreqdesc HbondVMfreqdesc HbondFCdesc Hbondindprintdesc];
  %Hbond_desc2=HbondVMintensitydesc;
  if numel(Hbond_desc1)>100 %256-6*26
    warning('rot09:toomuchcolumns','number of columns in Hbond datasection exceeds Excel limit! Decreased.')  
    Hbond_desc1=Hbond_desc1(1:100);
  end    
  xlswrite1spec(xlsfile,Hbond_desc1,'ribose','FA1');
  %xlswrite1spec(xlsfile,Hbond_desc2,'ribose',['FA' int2str(1+numel(Hbond_desc1))]);
  Hbond_data1=[nbonds' Hbonddata HbondLMfreq HbondVMfreq HbondFC Hbondindprint];
  if size(Hbond_data1,2)>100 %256-6*26
    warning('rot09:toomuchcolumns2','number of columns in Hbond datasection exceeds Excel limit! Decreased.')  
    Hbond_data1=Hbond_data1(:,1:100);
  end    
  %Hbond_data2=HbondVMintensity;
  xlswrite1spec(xlsfile,Hbond_data1,'ribose','FA3');
  %xlswrite1spec(xlsfile,Hbond_data2,'ribose',['FA' int2str(3+size(Hbond_data1,2))]);

if ~any(moltype==[232]) %zhr091227
  %symbols 'A'=65, 'F'=70
  calcHbondnum=cell(1,numel(Hbonddatadesc)+1);
  %for i=2:numel(Hbonddatadesc)+1
  for i=2:size(indHbond,1)
      calcHbondnum{i}=['=СЧЕТ(' indHbond(i,1:2) int2str(3) ':' indHbond(i,1:2) int2str(2+recnum) ')']; %{'=СЧЕТ(FA3:FA64)'}
  end
  xlswrite1spec(xlsfile,calcHbondnum,'ribose',['FA' int2str(recnum+3)])
end

  fdesc0=[{'hta2'},{'hta3'},{'hta4'},{'hta5'},{'hta0'}];
  fdesc=[{'tbeta'},{'tgamma'},{'tdelta'},{'tepsilon'},{'teta'},...
       {'tteta'},{'tchi'},{'tkappa2'},{'tnu0'},{'tnu1'},{'tnu2'},{'tnu3'},...
       {'tnu4'},{'sumofnu'},{'Pdeg'},{'numax'},cell(1,10),{'rglyc'},{'tau01'}];
  xlswrite1spec(xlsfile,[fdesc0 fdesc],'ribose',['AA' int2str(recnum+3)]);

  invoke(Excel.ActiveWorkbook,'Save');
  Excel.Quit
  Excel.delete
  clear Excel
  
end %endif flwritefile

cputime-time0   %#ok
toc
