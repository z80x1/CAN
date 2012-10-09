%rot28_2_importAIM_all: analyse AIM data in database for several molecules
%together
%
% Version 1.0    
% Last modified  R O Zhurakivsky 2010-12-18
% Created        R O Zhurakivsky 2010-12-18

tic

format compact

global pind
atomsind
pindsdef

clear AIM Hbondtype

%---------------------------------
%moltypes=[140 240 340 540]  %#ok
%theory='dftV2' %#ok

moltypes=[9 12 13 15 16]  %#ok
theory='dftV3' %#ok

%moltype=950    %#ok
usedpackage='Gaussian' %#ok
onlyoriginal=1;  % process db with only original conformations
version='01' %#ok

fl_load = 1 %#ok
flwritefile=1 %#ok
flwrite_bondlist=0 %#ok %output table with all CPs parameters
%---------------------------------

indir=[CD.datadir];
diaryfname0=[indir filesep 'logfile.rot28_2_all'];
diaryfname=diaryfname0;
for i=2:inf
  if ~(exist(diaryfname,'file')==2)
    break
  end
  diaryfname = [diaryfname0 int2str(i)];
end
diary(diaryfname)

ipC1=strcmpcellar(pind.labels,'pC1');
ipH12=strcmpcellar(pind.labels,'pH12');
ibO2=strcmpcellar(pind.labels,'bO2');
ibN3=strcmpcellar(pind.labels,'bN3');
ibC6=strcmpcellar(pind.labels,'bC6');
ibH6=strcmpcellar(pind.labels,'bH6');
ibC8=strcmpcellar(pind.labels,'bC8');
ibH8=strcmpcellar(pind.labels,'bH8');
ipO4=strcmpcellar(pind.labels,'pO4');
ipO2=strcmpcellar(pind.labels,'pO2');
ipH22=strcmpcellar(pind.labels,'pH22');
ipO5=strcmpcellar(pind.labels,'pO5');
ipH53=strcmpcellar(pind.labels,'pH53');

mols=[];
for n=1:numel(moltypes)
  mols=[mols '_r' int2str(moltypes(n))]; %#ok
end
xlsfile = mols;
if strcmp(usedpackage,'Gaussian')
  xlsfile = [xlsfile '_g'];
end
if ~strcmp(theory,'dft')
  xlsfile = [xlsfile '_' theory];
end
xlsfile=[xlsfile '_' version];
if onlyoriginal
    templ='_or';
    xlsfile = [xlsfile templ];
end
xlsfile=[CD.xlsdir filesep xlsfile '_AIMprops.xls'] %#ok

sdesc={};
hbondind=1;
allconfind=0; %сквозной индекс конформера, for sdesc filling

for n=1:numel(moltypes)
  moltype = moltypes(n);

  workdbname=['r' int2str(moltype)] %#ok

  if strcmp(usedpackage,'Gaussian')
    workdbname=[workdbname '_g'];
  end
  if ~strcmp(theory,'dft')
    workdbname=[workdbname '_' theory];
  end
  if onlyoriginal
      templ='_or';
      workdbname = [workdbname templ];
  end
  workdbname=[CD.dbdir filesep workdbname '.mat'] %#ok

  if fl_load
  	load(workdbname,'workdb')
  end
      
  recnum=numel(workdb);
  if n==1
    AIM.pinds=zeros(numel(moltypes)*recnum*4,3); %array size determination may be some better
    AIM.confind=zeros(numel(moltypes)*recnum*4,1);
  end
  
  for i=1:recnum
  	  allconfind = allconfind + 1;

      ms0 =  workdb(i);

      sdesc{allconfind} = ms0.prop.sdesc;
      workdb_all{allconfind} = ms0;

      if isfield(ms0,'AIM') && isfield(ms0.AIM,'pinds')
           AIMpind = ms0.AIM.pinds;
           if isfield(ms0.AIM,'nfrag')
               AIMnfrag = ms0.AIM.nfrag;
           else
               AIMnfrag = zeros(size(AIMpind));
           end
      else
           AIMpind=zeros(0,3);
           AIMnfrag = zeros(0,3);
      end
      hbondnumi=size(AIMpind,1);
      if hbondnumi
           AIM.pinds(hbondind:hbondind+hbondnumi-1,1:size(AIMpind,2)) = AIMpind;
           AIM.nfrag(hbondind:hbondind+hbondnumi-1,1:size(AIMpind,2)) = AIMnfrag;
           AIM.confind(hbondind:hbondind+hbondnumi-1) = allconfind;
           AIM.ro(hbondind:hbondind+hbondnumi-1) = ms0.AIM.ro;

           if isfield(ms0.AIM,'deltaroAIM2000')
               AIM.deltaro(hbondind:hbondind+hbondnumi-1) = -4*ms0.AIM.deltaroAIM2000;
           elseif isfield(ms0.AIM,'DelSqRho')
               AIM.deltaro(hbondind:hbondind+hbondnumi-1) = ms0.AIM.DelSqRho;
           end

           if isfield(ms0.AIM,'V')
               AIM.V(hbondind:hbondind+hbondnumi-1) = ms0.AIM.V;
           end
      end
           
      hbondind=hbondind+hbondnumi;

  %-------------------------------------------------------------
  end %i=1:recnum

end %n=1:numel(moltypes)
%-------------------------------------------------------------

if ~isfield(AIM,'ro')
    error('No ro field exists!')
elseif sum(isnan(AIM.ro))
    error('ro contains NaNs!')
end
if ~isfield(AIM,'deltaro')
    error('No deltaro field exists!')
elseif sum(isnan(AIM.deltaro))
    error('deltaro contains NaNs!')
end

AIM.pinds=AIM.pinds(1:hbondind-1,:);
AIM.confind=AIM.confind(1:hbondind-1);
AIM.Hbondtype=zeros(hbondind-1,1);

AIM.ABdist=zeros(1,hbondind-1);
AIM.HBdist=zeros(1,hbondind-1);
AIM.AHBang=zeros(1,hbondind-1);

%Hbondtype - array with H-bonds types data
uniquepinds=unique(AIM.pinds,'rows');
HBw2atoms=find(uniquepinds(:,3)==0);
if size(uniquepinds,2)==4
    HBw4atoms=find(uniquepinds(:,4)~=0);
else
    HBw4atoms=[];
end
HBw3atoms=setdiff(1:size(uniquepinds,1),[HBw2atoms; HBw4atoms]);
Hbondtype.pinds=[uniquepinds(HBw3atoms,:); uniquepinds(HBw4atoms,:); uniquepinds(HBw2atoms,:)];

for i=1:size(Hbondtype.pinds,1) %cycle by all Hbond's types

    curinds=find(sum(AIM.pinds==repmat(Hbondtype.pinds(i,:),hbondind-1,1),2)==size(AIM.pinds,2)); %indexes of Hbonds
    AIM.Hbondtype(curinds)=i; %type of Hbond

    conf_syn(i)=0;
    conf_anti(i)=0;
    conf_north(i)=0;
    conf_south(i)=0;
    conf_betagp(i)=0;
    conf_betat(i)=0;
    conf_betagm(i)=0;
    conf_gammagp(i)=0;
    conf_gammat(i)=0;
    conf_gammagm(i)=0;
    conf_deltagp(i)=0;
    conf_deltat(i)=0;
    conf_deltagm(i)=0;
    conf_epsilongp(i)=0;
    conf_epsilont(i)=0;
    conf_epsilongm(i)=0;
    if (moltype>100 && mod(moltype,100)==40) || moltype==7
        conf_etagp(i)=0;
        conf_etat(i)=0;
        conf_etagm(i)=0;
        conf_tetagp(i)=0;
        conf_tetat(i)=0;
        conf_tetagm(i)=0;
    end

    for j=1:numel(curinds)
      curconf=AIM.confind(curinds(j)); %number of current conformer

      ms0 = workdb_all{curconf};
      atom1ind=find(ms0.pind==Hbondtype.pinds(i,1));
      atom2ind=find(ms0.pind==Hbondtype.pinds(i,2));
      atom3ind=find(ms0.pind==Hbondtype.pinds(i,3));
      if isempty(atom1ind)
        error([ms0.prop.sdesc ': No atom ' pind.labels{Hbondtype.pinds(i,1)} ' exists in current molecule']);
      end
      if isempty(atom2ind)
        error([ms0.prop.sdesc ': No atom ' pind.labels{Hbondtype.pinds(i,2)} ' exists in current molecule']);
      end
%       if isempty(atom3ind)
%         error([ms0.prop.sdesc ': No atom ' pind.labels{Hbondtype.pinds(i,3)} ' exists in current molecule']);
%       end
          
      if ~isempty(atom3ind)
          AIM.ABdist(curinds(j))=adist(ms0,atom1ind,atom3ind);
          AIM.HBdist(curinds(j))=adist(ms0,atom2ind,atom3ind);
          AIM.AHBang(curinds(j))=valang(ms0,atom1ind,atom2ind,atom3ind);
      else
          AIM.ABdist(curinds(j))=0;
          AIM.HBdist(curinds(j))=adist(ms0,atom1ind,atom2ind);
          AIM.AHBang(curinds(j))=0;
      end

      if ( any(moltype==[7,8]) || any(moltype==[910,920,950]) )
      else
          ttt=ms0.prop.sdesc(end);
          if ttt=='S' || ttt=='T' || ttt=='V'
             conf_syn(i)=conf_syn(i)+1;
          elseif ttt=='A' || ttt=='B' || ttt=='C'
             conf_anti(i)=conf_anti(i)+1;
          else
             error(['Unknown conformation chintype detected: ' ms0.prop.sdesc]);
          end
      end
     
      if ms0.prop.Pdeg<90 || ms0.prop.Pdeg>270
         conf_north(i)=conf_north(i)+1;
      else
         conf_south(i)=conf_south(i)+1;
      end

      if ceil(ms0.prop.tbeta/120)==1  % (0,120]
         conf_betagp(i)=conf_betagp(i)+1;
      elseif ceil(ms0.prop.tbeta/120)==0 % (-120,0]
         conf_betagm(i)=conf_betagm(i)+1;
      else
         conf_betat(i)=conf_betat(i)+1;
      end

      if ceil(ms0.prop.tgamma/120)==1
         conf_gammagp(i)=conf_gammagp(i)+1;
      elseif ceil(ms0.prop.tgamma/120)==0
         conf_gammagm(i)=conf_gammagm(i)+1;
      else
         conf_gammat(i)=conf_gammat(i)+1;
      end

      if ceil(ms0.prop.tdelta/120)==1
         conf_deltagp(i)=conf_deltagp(i)+1;
      elseif ceil(ms0.prop.tdelta/120)==0
         conf_deltagm(i)=conf_deltagm(i)+1;
      else
         conf_deltat(i)=conf_deltat(i)+1;
      end

      if ceil(ms0.prop.tepsilon/120)==1
         conf_epsilongp(i)=conf_epsilongp(i)+1;
      elseif ceil(ms0.prop.tepsilon/120)==0
         conf_epsilongm(i)=conf_epsilongm(i)+1;
      else
         conf_epsilont(i)=conf_epsilont(i)+1;
      end

      if (moltype>100 && mod(moltype,100)==40) || moltype==7
          if ceil(ms0.prop.teta/120)==1
             conf_etagp(i)=conf_etagp(i)+1;
          elseif ceil(ms0.prop.teta/120)==0
             conf_etagm(i)=conf_etagm(i)+1;
          else
             conf_etat(i)=conf_etat(i)+1;
          end
          if ceil(ms0.prop.tteta/120)==1
             conf_tetagp(i)=conf_tetagp(i)+1;
          elseif ceil(ms0.prop.tteta/120)==0
             conf_tetagm(i)=conf_tetagm(i)+1;
          else
             conf_tetat(i)=conf_tetat(i)+1;
          end
      end

    end
      
    Hbondtype.num(i)=numel(curinds);

    Hbondtype.minABdist(i) = min(AIM.ABdist(curinds));
    Hbondtype.maxABdist(i) = max(AIM.ABdist(curinds));
    Hbondtype.avgABdist(i) = mean(AIM.ABdist(curinds));
    Hbondtype.stdABdist(i) = std(AIM.ABdist(curinds));

    Hbondtype.minHBdist(i) = min(AIM.HBdist(curinds));
    Hbondtype.maxHBdist(i) = max(AIM.HBdist(curinds));
    Hbondtype.avgHBdist(i) = mean(AIM.HBdist(curinds));
    Hbondtype.stdHBdist(i) = std(AIM.HBdist(curinds));

    Hbondtype.minAHBang(i) = min(AIM.AHBang(curinds));
    Hbondtype.maxAHBang(i) = max(AIM.AHBang(curinds));
    Hbondtype.avgAHBang(i) = mean(AIM.AHBang(curinds));
    Hbondtype.stdAHBang(i) = std(AIM.AHBang(curinds));

    Hbondtype.minro(i) = min(AIM.ro(curinds));
    Hbondtype.maxro(i) = max(AIM.ro(curinds));
    Hbondtype.avgro(i) = mean(AIM.ro(curinds));
    Hbondtype.stdro(i) = std(AIM.ro(curinds));

    Hbondtype.mindeltaro(i) = min(AIM.deltaro(curinds));
    Hbondtype.maxdeltaro(i) = max(AIM.deltaro(curinds));
    Hbondtype.avgdeltaro(i) = mean(AIM.deltaro(curinds));
    Hbondtype.stddeltaro(i) = std(AIM.deltaro(curinds));

    Hbondtype.minEhb(i) = min(abs(AIM.V(curinds))*0.5*CC.encoef); %minimal energy of H-bond, kcal/mol
    Hbondtype.maxEhb(i) = max(abs(AIM.V(curinds))*0.5*CC.encoef); %maximal energy of H-bond, kcal/mol

%    HbondVMfreq(curinds,Hbond_bH6_ind(curinds))

    if size(Hbondtype.pinds,2)==4 && Hbondtype.pinds(i,4)
        Hbondtype.title(i)={[pind.labels{Hbondtype.pinds(i,1)} pind.labels{Hbondtype.pinds(i,2)} '...' pind.labels{Hbondtype.pinds(i,3)} pind.labels{Hbondtype.pinds(i,4)}]};
    elseif Hbondtype.pinds(i,3)
        Hbondtype.title(i)={[pind.labels{Hbondtype.pinds(i,1)} pind.labels{Hbondtype.pinds(i,2)} '...' pind.labels{Hbondtype.pinds(i,3)}]};
    else
        Hbondtype.title(i)={[pind.labels{Hbondtype.pinds(i,1)}  '...' pind.labels{Hbondtype.pinds(i,2)}]};
    end
    
    AIM_type = Hbondtype.title{i}
    AIM_sdesc = sdesc(AIM.confind(curinds))
    AIM_HBdist = AIM.HBdist(curinds)

end %i=1:size(Hbondtype.pinds,1)

if any(AIM.pinds(find(AIM.ABdist==0),3)), error('Zero ABdist values detected.'), end %analyze only bonds with three nonzero pinds
if any(AIM.HBdist==0), error('Zero HBdist values detected.'), end
if any(AIM.pinds(find(AIM.ABdist==0),3)), error('Zero AHBang values detected.'), end

%-------------------------------------------------------
if flwritefile
  params.Bold = 1;
  params.Italic = 0;
  params.FontSize = 9;
  params.FontName = 'Verdana';
  params.ColorIndex = 1;
  params.NumberFormat = '@';

  try
    Excel = actxserver ('Excel.Application');
    File=fullfile(pwd,xlsfile);
    if ~exist(File,'file')
      ExcelWorkbook = Excel.workbooks.Add;
      ExcelWorkbook.SaveAs(File,1);
      ExcelWorkbook.Close(false);
    end
    invoke(Excel.Workbooks,'Open',File);

    Hbondtype2exportdesc =  [ {'H-bond name'}; {'number'};...
              {'minABdist'}; {'maxABdist'}; {'avgABdist'};{'stdABdist'}; {'stdABdist/avgABdist'};...
              {'minHBdist'}; {'maxHBdist'}; {'avgHBdist'};{'stdHBdist'}; {'stdHBdist/avgHBdist'};...
              {'minAHBang'}; {'maxAHBang'}; {'avgAHBang'};{'stdAHBang'}; {'stdAHBang/avgAHBang'};...
              {'minro'}; {'maxro'}; {'avgro'};{'stdro'}; {'stdro/avgro'};...
              {'mindeltaro'}; {'maxdeltaro'}; {'avgdeltaro'};{'stddeltaro'}; {'stddeltaro/avgdeltaro'};...
              {'minEhb, kcal/mol'}; {'maxEhb, kcal/mol'}; ...
                {'conf_syn'}; {'conf_anti'}; {'conf_north'}; {'conf_south'};...
                {'conf_betagp'}; {'conf_betat'}; {'conf_betagm'};...
                {'conf_gammagp'}; {'conf_gammat'}; {'conf_gammagm'};...
                {'conf_deltagp'}; {'conf_deltat'}; {'conf_deltagm'};...
                {'conf_epsilongp'}; {'conf_epsilont'}; {'conf_epsilongm'};
                ];
    if (moltype>100 && mod(moltype,100)==40) || moltype==7
        Hbondtype2exportdesc =  [ Hbondtype2exportdesc; {'conf_etagp'}; {'conf_etat'}; {'conf_etagm'};...
                                  {'conf_tetagp'}; {'conf_tetat'}; {'conf_tetagm'}];
    end
    xlswrite1spec(File,Hbondtype2exportdesc,'AIMavgdata','A1',params);

    l=1;
    xlswrite1spec(File,Hbondtype.title,'AIMavgdata',['B' int2str(l)],params); l=l+1;
    params.Bold = 0;
    params.Italic = 0;
    params.NumberFormat = '0';
    xlswrite1spec(File,Hbondtype.num,'AIMavgdata',['B' int2str(l)],params); l=l+1;
    params.NumberFormat = '0.000';
    xlswrite1spec(File,Hbondtype.minABdist,'AIMavgdata',['B' int2str(l)],params); l=l+1;
    xlswrite1spec(File,Hbondtype.maxABdist,'AIMavgdata',['B' int2str(l)],params); l=l+1;
    xlswrite1spec(File,Hbondtype.avgABdist,'AIMavgdata',['B' int2str(l)],params); l=l+1;
    xlswrite1spec(File,Hbondtype.stdABdist,'AIMavgdata',['B' int2str(l)],params); l=l+1;
    xlswrite1spec(File,Hbondtype.stdABdist./Hbondtype.avgABdist,'AIMavgdata',['B' int2str(l)],params); l=l+1;
    xlswrite1spec(File,Hbondtype.minHBdist,'AIMavgdata',['B' int2str(l)],params); l=l+1;
    xlswrite1spec(File,Hbondtype.maxHBdist,'AIMavgdata',['B' int2str(l)],params); l=l+1;
    xlswrite1spec(File,Hbondtype.avgHBdist,'AIMavgdata',['B' int2str(l)],params); l=l+1;
    xlswrite1spec(File,Hbondtype.stdHBdist,'AIMavgdata',['B' int2str(l)],params); l=l+1;
    xlswrite1spec(File,Hbondtype.stdHBdist./Hbondtype.avgHBdist,'AIMavgdata',['B' int2str(l)],params); l=l+1;
    params.NumberFormat = '0.0';
    xlswrite1spec(File,Hbondtype.minAHBang,'AIMavgdata',['B' int2str(l)],params); l=l+1;
    xlswrite1spec(File,Hbondtype.maxAHBang,'AIMavgdata',['B' int2str(l)],params); l=l+1;
    xlswrite1spec(File,Hbondtype.avgAHBang,'AIMavgdata',['B' int2str(l)],params); l=l+1;
    xlswrite1spec(File,Hbondtype.stdAHBang,'AIMavgdata',['B' int2str(l)],params); l=l+1;
    xlswrite1spec(File,Hbondtype.stdAHBang./Hbondtype.avgAHBang,'AIMavgdata',['B' int2str(l)],params); l=l+1;
    params.NumberFormat = '0.000';
    xlswrite1spec(File,Hbondtype.minro,'AIMavgdata',['B' int2str(l)],params); l=l+1;
    xlswrite1spec(File,Hbondtype.maxro,'AIMavgdata',['B' int2str(l)],params); l=l+1;
    xlswrite1spec(File,Hbondtype.avgro,'AIMavgdata',['B' int2str(l)],params); l=l+1;
    xlswrite1spec(File,Hbondtype.stdro,'AIMavgdata',['B' int2str(l)],params); l=l+1;
    xlswrite1spec(File,Hbondtype.stdro./Hbondtype.avgro,'AIMavgdata',['B' int2str(l)],params); l=l+1;
    xlswrite1spec(File,Hbondtype.mindeltaro,'AIMavgdata',['B' int2str(l)],params); l=l+1;
    xlswrite1spec(File,Hbondtype.maxdeltaro,'AIMavgdata',['B' int2str(l)],params); l=l+1;
    xlswrite1spec(File,Hbondtype.avgdeltaro,'AIMavgdata',['B' int2str(l)],params); l=l+1;
    xlswrite1spec(File,Hbondtype.stddeltaro,'AIMavgdata',['B' int2str(l)],params); l=l+1;
    xlswrite1spec(File,Hbondtype.stddeltaro./Hbondtype.avgdeltaro,'AIMavgdata',['B' int2str(l)],params); l=l+1;

    params.NumberFormat = '0.0';
    xlswrite1spec(File,Hbondtype.minEhb,'AIMavgdata',['B' int2str(l)],params); l=l+1;
    xlswrite1spec(File,Hbondtype.maxEhb,'AIMavgdata',['B' int2str(l)],params); l=l+1;

    params.NumberFormat = '0';
    confdata = [conf_syn; conf_anti; conf_north; conf_south;...
        conf_betagp; conf_betat; conf_betagm;...
        conf_gammagp; conf_gammat; conf_gammagm;...
        conf_deltagp; conf_deltat; conf_deltagm;...
        conf_epsilongp; conf_epsilont; conf_epsilongm];
    if (moltype>100 && mod(moltype,100)==40) || moltype==7
        confdata = [confdata;  conf_etagp; conf_etat; conf_etagm;...
                    conf_tetagp; conf_tetat; conf_tetagm];
    end
    confdata = num2cell(confdata);
    xlswrite1spec(File,confdata,'AIMavgdata',['B' int2str(l)],params); l=l+size(confdata,1);


    %---ro, deltaro & Hbond's geometry parameters correlations

    inds=find(AIM.pinds(:,3)~=0); %exclude bonds with 2 atoms (VdV interactions) from correlation analysis
    ro=AIM.ro(inds);
    deltaro=AIM.deltaro(inds);
    ABdist=AIM.ABdist(inds);
    HBdist=AIM.HBdist(inds);
    AHBang=AIM.AHBang(inds);
    
    corrdata=[ro' deltaro' ABdist' HBdist' AHBang'];
    corrdatadesc=[{'ro'}; {'deltaro'}; {'ABdist'}; {'HBdist'}; {'AHBang'}];
    [r,p] = corrcoef(corrdata);
    %r(find(p>0.000001))=NaN;
    %ii=find(abs(r)<0.7);
    rr=num2cell(r);
    %rr(ii)={'-'};
    %rr=reshape(rr,size(p));

    ttt=repmat({''},1,size(rr,1));
    tcorrcoefmatr=[corrdatadesc ttt' rr];
    tcorrcoefmatr=[{''} {''} corrdatadesc'; {''} {''} ttt; tcorrcoefmatr];
    %----------------------------------------------

    params.NumberFormat = '0.00';
    xlswrite1spec(File,tcorrcoefmatr,'AIMavgdata',['A' int2str(l+1)],params); l=l+size(tcorrcoefmatr,1);


%    xlswrite1spec(File,'Data for all Hbonds detected by AIM','AIMavgdata',['A' int2str(l+1)],params);
%    xlswrite1spec(File,corrdatadesc,'AIMavgdata',['A' int2str(l+2)],params);
%    xlswrite1spec(File,corrdata','AIMavgdata',['B' int2str(l+2)],params);
%    %maybe too many columns for Excel
    

    invoke(Excel.ActiveWorkbook,'Save');
    Excel.Quit
    Excel.delete
    clear Excel
    
  catch
    Excel.Quit
    Excel.delete
    clear Excel
  end

end %flwritefile



toc
diary off
clear workdb_all