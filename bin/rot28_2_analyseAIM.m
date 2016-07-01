%rot28_2_importAIM: analyse AIM data in database
%
% Version 1.0    
% Last modified  R O Zhurakivsky 2009-04-25
% Created        R O Zhurakivsky 2006-?-?

%070829 added analisys of population of conformers with each type of Hbonds

tic
clear 
format compact

global pind
atomsind
pindsdef

%---------------------------------
T=298.15        %#ok
moltype=9  %#ok
usedpackage='Gaussian' %#ok
theory='dftV3' %#ok
onlyoriginal=1;  % process db with only original conformations
version='01' %#ok
flwritefile=1 %#ok

flanalyzeC1H12=0;
flanalyzeC6H6=0;
flanalyzeC8H8=0;
flanalyze_pO2HbO2=0;
flanalyze_pO5HbO2=0;
flanalyze_pO2HbN3=0;
flanalyze_pO5HbN3=0;
%---------------------------------

indir=[CD.datadir filesep 'r' int2str(moltype)];
diaryfname0=[indir filesep 'logfile.rot28_2'];
diaryfname=diaryfname0;
for i=2:inf
  if ~(exist(diaryfname,'file')==2)
    break
  end
  diaryfname = [diaryfname0 int2str(i)];
end
diary(diaryfname)

workdbname=['r' int2str(moltype)] %#ok
xlsfile = workdbname;

if strcmp(usedpackage,'Gaussian')
  workdbname=[workdbname '_g'];
  xlsfile = [xlsfile '_g'];
end
if ~strcmp(theory,'dft')
  workdbname=[workdbname '_' theory];
  xlsfile = [xlsfile '_' theory];
end
xlsfile=[xlsfile '_' version];
if onlyoriginal
    templ='_or';
    workdbname = [workdbname templ];
    xlsfile = [xlsfile templ];
end
workdbname=[CD.dbdir filesep workdbname '.mat'] %#ok
xlsfile=[CD.xlsdir filesep xlsfile '.xls'] %#ok


load(workdbname,'workdb')

recnum=numel(workdb);

sdesc={};
AIM.pinds=zeros(recnum*4,3);
AIM.confind=zeros(recnum*4,1);
hbondind=1;
conf_wo_O5_=0; %number of conformations without O5' atom in Hbonds

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

conf_onlygeomHbondsexists=0;
pC1pH12bO2notAIManti.sdesc={};
pC1pH12bO2notAIManti.AHdist=[];
pC1pH12bO2notAIManti.ABdist=[];
pC1pH12bO2notAIManti.HBdist=[];
pC1pH12bO2notAIManti.AHBang=[];
bC6bH6pO4notAIManti.sdesc={};
bC6bH6pO4notAIManti.AHdist=[];
bC6bH6pO4notAIManti.ABdist=[];
bC6bH6pO4notAIManti.HBdist=[];
bC6bH6pO4notAIManti.AHBang=[];

bC6bH6pO4pO5notAIManti.sdesc={};
bC6bH6pO4pO5notAIManti.AHdist=[];
bC6bH6pO4pO5notAIManti.VMfreq=[];
bC6bH6pO4pO5notAIManti.VMint=[];
bC6bH6pO4_AIManti.sdesc={};
bC6bH6pO4_AIManti.AHdist=[];
bC6bH6pO4_AIManti.VMfreq=[];
bC6bH6pO4_AIManti.VMint=[];
bC6bH6pO5_AIManti.sdesc={};
bC6bH6pO5_AIManti.AHdist=[];
bC6bH6pO5_AIManti.VMfreq=[];
bC6bH6pO5_AIManti.VMint=[];
bC6bH6pO4pO5_AIManti.sdesc={};
bC6bH6pO4pO5_AIManti.AHdist=[];
bC6bH6pO4pO5_AIManti.VMfreq=[];
bC6bH6pO4pO5_AIManti.VMint=[];

pO2pH22bO2notAIManti.sdesc={};
pO2pH22bO2notAIManti.AHdist=[];
pO2pH22bO2notAIManti.VMfreq=[];
pO2pH22bO2notAIManti.VMint=[];
pO2pH22bO2_AIManti.sdesc={};
pO2pH22bO2_AIManti.AHdist=[];
pO2pH22bO2_AIManti.VMfreq=[];
pO2pH22bO2_AIManti.VMint=[];

pO5pH53bO2notAIMsyn.sdesc={};
pO5pH53bO2notAIMsyn.AHdist=[];
pO5pH53bO2notAIMsyn.VMfreq=[];
pO5pH53bO2notAIMsyn.VMint=[];
pO5pH53bO2_AIMsyn.sdesc={};
pO5pH53bO2_AIMsyn.AHdist=[];
pO5pH53bO2_AIMsyn.VMfreq=[];
pO5pH53bO2_AIMsyn.VMint=[];

pO2pH22bN3notAIManti.sdesc={};
pO2pH22bN3notAIManti.AHdist=[];
pO2pH22bN3notAIManti.VMfreq=[];
pO2pH22bN3notAIManti.VMint=[];
pO2pH22bN3_AIManti.sdesc={};
pO2pH22bN3_AIManti.AHdist=[];
pO2pH22bN3_AIManti.VMfreq=[];
pO2pH22bN3_AIManti.VMint=[];

pO5pH53bN3notAIMsyn.sdesc={};
pO5pH53bN3notAIMsyn.AHdist=[];
pO5pH53bN3notAIMsyn.VMfreq=[];
pO5pH53bN3notAIMsyn.VMint=[];
pO5pH53bN3_AIMsyn.sdesc={};
pO5pH53bN3_AIMsyn.AHdist=[];
pO5pH53bN3_AIMsyn.VMfreq=[];
pO5pH53bN3_AIMsyn.VMint=[];

bC8bH8pO5notAIManti.sdesc={};
bC8bH8pO5notAIManti.AHdist=[];
bC8bH8pO5notAIManti.VMfreq=[];
bC8bH8pO5notAIManti.VMint=[];
bC8bH8pO5_AIManti.sdesc={};
bC8bH8pO5_AIManti.AHdist=[];
bC8bH8pO5_AIManti.VMfreq=[];
bC8bH8pO5_AIManti.VMint=[];

modeind=[];

if isfield(workdb(1).gaussian,'MP2_311__G2dfpd')
    energyfield='MP2_311__G2dfpd';
elseif isfield(workdb(1).gaussian,'MP2_6311__G2dfpd')
    energyfield='MP2_6311__G2dfpd';
elseif isfield(workdb(1).gaussian,'MP2_6311__Gdp')
    energyfield='MP2_6311__Gdp';
else
    error('energy field not found');
end
disp(['energy field ' energyfield ' detected'])
energy=[];

for i=1:recnum

    ms0 =  workdb(i);
    sdesc{i} = ms0.prop.sdesc;

%     if ~numel(ms0.AIM.desc)
%         continue
%     end
    if isfield(ms0.gaussian,'T')
       tind = find(ms0.gaussian.T==T);
       if isempty(tind), error(['Desired temperature T=' num2str(T) 'K is not found']), end  
    else %???
       tind=1;
    end
    if isfield(ms0.gaussian,energyfield)
         energy(end+1) = ms0.gaussian.(energyfield) + ms0.gaussian.GEC(tind)/CC.encoef;
    else
         error(['energy field ' energyfield ' in record # ' int2str(i) ' is not found']);
    end

    if ~isempty(ms0.AIM) && isfield(ms0.AIM,'pinds')
         AIMpind = ms0.AIM.pinds;
    else
         AIMpind=zeros(0,3);
    end
    hbondnumi=size(AIMpind,1);
    if hbondnumi
         AIM.pinds(hbondind:hbondind+hbondnumi-1,1:size(AIMpind,2)) = AIMpind;
         AIM.confind(hbondind:hbondind+hbondnumi-1) = i;
         AIM.ro(hbondind:hbondind+hbondnumi-1) = ms0.AIM.ro;

         if isfield(ms0.AIM,'deltaroAIM2000')
             AIM.deltaro(hbondind:hbondind+hbondnumi-1) = -4*ms0.AIM.deltaroAIM2000;
         elseif isfield(ms0.AIM,'DelSqRho')
             AIM.deltaro(hbondind:hbondind+hbondnumi-1) = ms0.AIM.DelSqRho;
         end
    end
         
    if any(any(AIMpind==ipO5)) %???
       conf_wo_O5_ = conf_wo_O5_ + 1;
    end

    if 0
         geompind = [ms0.pind(ms0.prop.HbondO1ind) ms0.pind(ms0.prop.HbondHind) ms0.pind(ms0.prop.HbondO2ind)];
         
         [onlygeomHbonds,II]=setdiff(geompind,AIMpind,'rows');

         if ~isempty(onlygeomHbonds)
             conf_onlygeomHbondsexists=conf_onlygeomHbondsexists+1;
             AHBang=0; ABdist=0; HBdist=0;
             for jj=1:numel(II)
               AHBang(jj)=valang(ms0,ms0.prop.HbondO1ind(II(jj)),ms0.prop.HbondHind(II(jj)),ms0.prop.HbondO2ind(II(jj)));
               ABdist(jj)=adist(ms0,ms0.prop.HbondO1ind(II(jj)),ms0.prop.HbondO2ind(II(jj)));
               HBdist(jj)=adist(ms0,ms0.prop.HbondHind(II(jj)),ms0.prop.HbondO2ind(II(jj)));
             end 
             disp(['#' int2str(conf_onlygeomHbondsexists) ' ' ms0.prop.sdesc ': AHBang ' num2str(AHBang,'%4.1f ') ', ABdist ' num2str(ABdist,'%4.2f ') ', HBdist ' num2str(HBdist,'%4.2f ')])
             disp(reshape(pind.labels(onlygeomHbonds),size(onlygeomHbonds)))
             disp(' ')
         end
    end
     
    hbondind=hbondind+hbondnumi;

%-------------------------------------------------------------
    if flanalyzeC1H12
     if ~any(all(repmat([ipC1 ipH12 ibO2]',1,size(AIMpind,1))==AIMpind')) && ms0.prop.sdesc(end)=='A'
     %this is anti conformation  without C1'H...O2 bond detected by AIM
        atom1ind=find(ms0.pind==ipC1);
        atom2ind=find(ms0.pind==ipH12);
        atom3ind=find(ms0.pind==ibO2);
        pC1pH12bO2notAIManti.ABdist(end+1)=adist(ms0,atom1ind,atom3ind);
        pC1pH12bO2notAIManti.HBdist(end+1)=adist(ms0,atom2ind,atom3ind);
        pC1pH12bO2notAIManti.AHBang(end+1)=valang(ms0,atom1ind,atom2ind,atom3ind);
     end
    end
    
    if flanalyzeC6H6
        if ~any(all(repmat([ibC6 ibH6 ipO4]',1,size(AIMpind,1))==AIMpind')) && ...
             ~any(all(repmat([ibC6 ibH6 ipO5]',1,size(AIMpind,1))==AIMpind')) && ms0.prop.sdesc(end)=='A'
          %this is anti conformation  without C6H...O4' & C6H...O5' bonds detected by AIM
             atom1ind=find(ms0.pind==ibC6);
             atom2ind=find(ms0.pind==ibH6);
        %     if isempty(atom1ind) || isempty(atom2ind), error('empty indexes detected'), end;
             bC6bH6pO4pO5notAIManti.sdesc(end+1)={ms0.prop.sdesc};
             bC6bH6pO4pO5notAIManti.AHdist(end+1)=adist(ms0,atom1ind,atom2ind);
             modeind(end+1)=ms0.prop.HbondVMind(ms0.pind(strcmpcellar(ms0.labels,'H'))==ibH6);
             bC6bH6pO4pO5notAIManti.VMfreq(end+1)=ms0.freq.freq(modeind(end));
             bC6bH6pO4pO5notAIManti.VMint(end+1)=ms0.freq.intencity(modeind(end));
        end

        if any(all(repmat([ibC6 ibH6 ipO4]',1,size(AIMpind,1))==AIMpind')) && ...
             ~any(all(repmat([ibC6 ibH6 ipO5]',1,size(AIMpind,1))==AIMpind')) && ms0.prop.sdesc(end)=='A'
          %this is anti conformation  with only C6H...O4' bond detected by AIM
             atom1ind=find(ms0.pind==ibC6);
             atom2ind=find(ms0.pind==ibH6);
       %      if isempty(atom1ind) || isempty(atom2ind), error('empty indexes detected'), end;
             bC6bH6pO4_AIManti.sdesc(end+1)={ms0.prop.sdesc};
             bC6bH6pO4_AIManti.AHdist(end+1)=adist(ms0,atom1ind,atom2ind);
             modeind(end+1)=ms0.prop.HbondVMind(ms0.pind(strcmpcellar(ms0.labels,'H'))==ibH6);
             bC6bH6pO4_AIManti.VMfreq(end+1)=ms0.freq.freq(modeind(end));
             bC6bH6pO4_AIManti.VMint(end+1)=ms0.freq.intencity(modeind(end));
        end
        if ~any(all(repmat([ibC6 ibH6 ipO4]',1,size(AIMpind,1))==AIMpind')) && ...
             any(all(repmat([ibC6 ibH6 ipO5]',1,size(AIMpind,1))==AIMpind')) && ms0.prop.sdesc(end)=='A'
          %this is anti conformation  with only C6H...O5' bond detected by AIM
             atom1ind=find(ms0.pind==ibC6);
             atom2ind=find(ms0.pind==ibH6);
      %       if isempty(atom1ind) || isempty(atom2ind), error('empty indexes detected'), end;
             bC6bH6pO5_AIManti.sdesc(end+1)={ms0.prop.sdesc};
             bC6bH6pO5_AIManti.AHdist(end+1)=adist(ms0,atom1ind,atom2ind);
             modeind(end+1)=ms0.prop.HbondVMind(ms0.pind(strcmpcellar(ms0.labels,'H'))==ibH6);
             bC6bH6pO5_AIManti.VMfreq(end+1)=ms0.freq.freq(modeind(end));
             bC6bH6pO5_AIManti.VMint(end+1)=ms0.freq.intencity(modeind(end));
        end
        if any(all(repmat([ibC6 ibH6 ipO4]',1,size(AIMpind,1))==AIMpind')) && ...
             any(all(repmat([ibC6 ibH6 ipO5]',1,size(AIMpind,1))==AIMpind')) && ms0.prop.sdesc(end)=='A'
          %this is anti conformation  with C6H...O4' & C6H...O5' bonds detected by AIM
             atom1ind=find(ms0.pind==ibC6);
             atom2ind=find(ms0.pind==ibH6);
     %        if isempty(atom1ind) || isempty(atom2ind), error('empty indexes detected'), end;
             bC6bH6pO4pO5_AIManti.sdesc(end+1)={ms0.prop.sdesc};
             bC6bH6pO4pO5_AIManti.AHdist(end+1)=adist(ms0,atom1ind,atom2ind);
             modeind(end+1)=ms0.prop.HbondVMind(ms0.pind(strcmpcellar(ms0.labels,'H'))==ibH6);
             bC6bH6pO4pO5_AIManti.VMfreq(end+1)=ms0.freq.freq(modeind(end));
             bC6bH6pO4pO5_AIManti.VMint(end+1)=ms0.freq.intencity(modeind(end));
        end
    end %flanalyzeC6H6

    if flanalyze_pO2HbO2
        if ~any(all(repmat([ipO2 ipH22 ibO2]',1,size(AIMpind,1))==AIMpind')) && ms0.prop.sdesc(end)=='A'
          %this is anti conformation  without O2'H...O2 bonds detected by AIM
             atom1ind=find(ms0.pind==ipO2);
             atom2ind=find(ms0.pind==ipH22);
    %         if isempty(atom1ind) || isempty(atom2ind), error('empty indexes detected'), end;
             pO2pH22bO2notAIManti.sdesc(end+1)={ms0.prop.sdesc};
             pO2pH22bO2notAIManti.AHdist(end+1)=adist(ms0,atom1ind,atom2ind);
             modeind(end+1)=ms0.prop.HbondVMind(ms0.pind(strcmpcellar(ms0.labels,'H'))==ipH22);
             pO2pH22bO2notAIManti.VMfreq(end+1)=ms0.freq.freq(modeind(end));
             pO2pH22bO2notAIManti.VMint(end+1)=ms0.freq.intencity(modeind(end));
        end

        if any(all(repmat([ipO2 ipH22 ibO2]',1,size(AIMpind,1))==AIMpind')) && ms0.prop.sdesc(end)=='A'
          %this is anti conformation  with only O2'H...O2 bond detected by AIM
             atom1ind=find(ms0.pind==ipO2);
             atom2ind=find(ms0.pind==ipH22);
     %        if isempty(atom1ind) || isempty(atom2ind), error('empty indexes detected'), end;
             pO2pH22bO2_AIManti.sdesc(end+1)={ms0.prop.sdesc};
             pO2pH22bO2_AIManti.AHdist(end+1)=adist(ms0,atom1ind,atom2ind);
             modeind(end+1)=ms0.prop.HbondVMind(ms0.pind(strcmpcellar(ms0.labels,'H'))==ipH22);
             pO2pH22bO2_AIManti.VMfreq(end+1)=ms0.freq.freq(modeind(end));
             pO2pH22bO2_AIManti.VMint(end+1)=ms0.freq.intencity(modeind(end));
        end
    end %flanalyze_pO2HbO2
    if flanalyze_pO5HbO2
      if ~any(all(repmat([ipO5 ipH53 ibO2]',1,size(AIMpind,1))==AIMpind')) && ms0.prop.sdesc(end)=='S'
      %this is syn conformation  without O5'H...O2 bonds detected by AIM
         atom1ind=find(ms0.pind==ipO5);
         atom2ind=find(ms0.pind==ipH53);
%         if isempty(atom1ind) || isempty(atom2ind), error('empty indexes detected'), end;
         pO5pH53bO2notAIMsyn.sdesc(end+1)={ms0.prop.sdesc};
         pO5pH53bO2notAIMsyn.AHdist(end+1)=adist(ms0,atom1ind,atom2ind);
         modeind(end+1)=ms0.prop.HbondVMind(ms0.pind(strcmpcellar(ms0.labels,'H'))==ipH53);
         pO5pH53bO2notAIMsyn.VMfreq(end+1)=ms0.freq.freq(modeind(end));
         pO5pH53bO2notAIMsyn.VMint(end+1)=ms0.freq.intencity(modeind(end));
      end

      if any(all(repmat([ipO5 ipH53 ibO2]',1,size(AIMpind,1))==AIMpind')) && ms0.prop.sdesc(end)=='S'
      %this is syn conformation  with only O5'H...O2 bond detected by AIM
         atom1ind=find(ms0.pind==ipO5);
         atom2ind=find(ms0.pind==ipH53);
%         if isempty(atom1ind) || isempty(atom2ind), error('empty indexes detected'), end;
         pO5pH53bO2_AIMsyn.sdesc(end+1)={ms0.prop.sdesc};
         pO5pH53bO2_AIMsyn.AHdist(end+1)=adist(ms0,atom1ind,atom2ind);
         modeind(end+1)=ms0.prop.HbondVMind(ms0.pind(strcmpcellar(ms0.labels,'H'))==ipH53);
         pO5pH53bO2_AIMsyn.VMfreq(end+1)=ms0.freq.freq(modeind(end));
         pO5pH53bO2_AIMsyn.VMint(end+1)=ms0.freq.intencity(modeind(end));
      end
    end %flanalyze_pO5HbO2

    if flanalyze_pO2HbN3
      if ~any(all(repmat([ipO2 ipH22 ibN3]',1,size(AIMpind,1))==AIMpind')) && ms0.prop.sdesc(end)=='A'
      %this is anti conformation  without O2'H...N3 bonds detected by AIM
         atom1ind=find(ms0.pind==ipO2);
         atom2ind=find(ms0.pind==ipH22);
%         if isempty(atom1ind) || isempty(atom2ind), error('empty indexes detected'), end;
         pO2pH22bN3notAIManti.sdesc(end+1)={ms0.prop.sdesc};
         pO2pH22bN3notAIManti.AHdist(end+1)=adist(ms0,atom1ind,atom2ind);
         modeind(end+1)=ms0.prop.HbondVMind(ms0.pind(strcmpcellar(ms0.labels,'H'))==ipH22);
         pO2pH22bN3notAIManti.VMfreq(end+1)=ms0.freq.freq(modeind(end));
         pO2pH22bN3notAIManti.VMint(end+1)=ms0.freq.intencity(modeind(end));
      end

      if any(all(repmat([ipO2 ipH22 ibN3]',1,size(AIMpind,1))==AIMpind')) && ms0.prop.sdesc(end)=='A'
      %this is anti conformation  with only O2'H...N3 bond detected by AIM
         atom1ind=find(ms0.pind==ipO2);
         atom2ind=find(ms0.pind==ipH22);
%         if isempty(atom1ind) || isempty(atom2ind), error('empty indexes detected'), end;
         pO2pH22bN3_AIManti.sdesc(end+1)={ms0.prop.sdesc};
         pO2pH22bN3_AIManti.AHdist(end+1)=adist(ms0,atom1ind,atom2ind);
         modeind(end+1)=ms0.prop.HbondVMind(ms0.pind(strcmpcellar(ms0.labels,'H'))==ipH22);
         pO2pH22bN3_AIManti.VMfreq(end+1)=ms0.freq.freq(modeind(end));
         pO2pH22bN3_AIManti.VMint(end+1)=ms0.freq.intencity(modeind(end));
      end
    end %flanalyze_pO2HbN3
    if flanalyze_pO5HbN3
      if ~any(all(repmat([ipO5 ipH53 ibN3]',1,size(AIMpind,1))==AIMpind')) && ms0.prop.sdesc(end)=='S'
      %this is syn conformation  without O5'H...N3 bonds detected by AIM
         atom1ind=find(ms0.pind==ipO5);
         atom2ind=find(ms0.pind==ipH53);
%         if isempty(atom1ind) || isempty(atom2ind), error('empty indexes detected'), end;
         pO5pH53bN3notAIMsyn.sdesc(end+1)={ms0.prop.sdesc};
         pO5pH53bN3notAIMsyn.AHdist(end+1)=adist(ms0,atom1ind,atom2ind);
         modeind(end+1)=ms0.prop.HbondVMind(ms0.pind(strcmpcellar(ms0.labels,'H'))==ipH53);
         pO5pH53bN3notAIMsyn.VMfreq(end+1)=ms0.freq.freq(modeind(end));
         pO5pH53bN3notAIMsyn.VMint(end+1)=ms0.freq.intencity(modeind(end));
      end

      if any(all(repmat([ipO5 ipH53 ibN3]',1,size(AIMpind,1))==AIMpind')) && ms0.prop.sdesc(end)=='S'
      %this is syn conformation  with only O5'H...N3 bond detected by AIM
         atom1ind=find(ms0.pind==ipO5);
         atom2ind=find(ms0.pind==ipH53);
 %        if isempty(atom1ind) || isempty(atom2ind), error('empty indexes detected'), end;
         pO5pH53bN3_AIMsyn.sdesc(end+1)={ms0.prop.sdesc};
         pO5pH53bN3_AIMsyn.AHdist(end+1)=adist(ms0,atom1ind,atom2ind);
         modeind(end+1)=ms0.prop.HbondVMind(ms0.pind(strcmpcellar(ms0.labels,'H'))==ipH53);
         pO5pH53bN3_AIMsyn.VMfreq(end+1)=ms0.freq.freq(modeind(end));
         pO5pH53bN3_AIMsyn.VMint(end+1)=ms0.freq.intencity(modeind(end));
      end
    end %flanalyze_pO5HbN3
%     HbondVMfreq(i,:) = ms0.freq.freq(ms0.prop.HbondVMind);
%     HbondVMintensity(i,:) = ms0.freq.intensity(ms0.prop.HbondVMind);
%     Hbond_bH6_ind(i) = find(ms0.pinds(strcellcmp(ms0.labels,'H')==ibH6));
%    HbondVMfreq(curinds,Hbond_bH6_ind(curinds))

    if flanalyzeC8H8
     if ~any(all(repmat([ibC8 ibH8 ipO5]',1,size(AIMpind,1))==AIMpind')) && ms0.prop.sdesc(end)=='A'
     %this is anti conformation  without C8H...O5' bonds detected by AIM
        atom1ind=find(ms0.pind==ibC8);
        atom2ind=find(ms0.pind==ibH8);
%        if isempty(atom1ind) || isempty(atom2ind), error('empty indexes detected'), end;
        bC8bH8pO5notAIManti.sdesc(end+1)={ms0.prop.sdesc};
        bC8bH8pO5notAIManti.AHdist(end+1)=adist(ms0,atom1ind,atom2ind);
        modeind(end+1)=ms0.prop.HbondVMind(ms0.pind(strcmpcellar(ms0.labels,'H'))==ibH8);
        bC8bH8pO5notAIManti.VMfreq(end+1)=ms0.freq.freq(modeind(end));
        bC8bH8pO5notAIManti.VMint(end+1)=ms0.freq.intencity(modeind(end));
     end

     if any(all(repmat([ibC8 ibH8 ipO5]',1,size(AIMpind,1))==AIMpind')) && ms0.prop.sdesc(end)=='A'
     %this is anti conformation  with only C8H...O5' bond detected by AIM
        atom1ind=find(ms0.pind==ibC8);
        atom2ind=find(ms0.pind==ibH8);
%        if isempty(atom1ind) || isempty(atom2ind), error('empty indexes detected'), end;
        bC8bH8pO5_AIManti.sdesc(end+1)={ms0.prop.sdesc};
        bC8bH8pO5_AIManti.AHdist(end+1)=adist(ms0,atom1ind,atom2ind);
        modeind(end+1)=ms0.prop.HbondVMind(ms0.pind(strcmpcellar(ms0.labels,'H'))==ibH8);
        bC8bH8pO5_AIManti.VMfreq(end+1)=ms0.freq.freq(modeind(end));
        bC8bH8pO5_AIManti.VMint(end+1)=ms0.freq.intencity(modeind(end));
     end
    end%flanalyzeC8H8
end
%-------------------------------------------------------------

energy=(energy-min(energy))*CC.hartree;
normcoef=1/(sum(exp(-energy/CC.k/T)));
population=normcoef*exp(-energy/CC.k/T); %#ok

if sum(isnan(AIM.ro))
    error('ro contains NaNs!')
end
if sum(isnan(AIM.deltaro))
    error('deltaro contains NaNs!')
end

iplotsigns=0;
plotsigns='.ox+*sdv^<>ph';

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
    conf_Htype_popul(i)=0;%population of all conformers that have this type of Hbond

    for j=1:numel(curinds)
      curconf=AIM.confind(curinds(j)); %number of current conformer
      atom1ind=find(workdb(curconf).pind==Hbondtype.pinds(i,1));
      atom2ind=find(workdb(curconf).pind==Hbondtype.pinds(i,2));
      atom3ind=find(workdb(curconf).pind==Hbondtype.pinds(i,3));
      if isempty(atom1ind)
        error([workdb(curconf).prop.sdesc ': No atom ' pind.labels{Hbondtype.pinds(i,1)} ' exists in current molecule']);
      end
      if isempty(atom2ind)
        error([workdb(curconf).prop.sdesc ': No atom ' pind.labels{Hbondtype.pinds(i,2)} ' exists in current molecule']);
      end
%       if isempty(atom3ind)
%         error([workdb(curconf).prop.sdesc ': No atom ' pind.labels{Hbondtype.pinds(i,3)} ' exists in current molecule']);
%       end
          
      if ~isempty(atom3ind)
          AIM.ABdist(curinds(j))=adist(workdb(curconf),atom1ind(1),atom3ind(1));
          AIM.HBdist(curinds(j))=adist(workdb(curconf),atom2ind(1),atom3ind(1));
          AIM.AHBang(curinds(j))=valang(workdb(curconf),atom1ind(1),atom2ind(1),atom3ind(1));
      else
          AIM.ABdist(curinds(j))=0;
          AIM.HBdist(curinds(j))=adist(workdb(curconf),atom1ind,atom2ind);
          AIM.AHBang(curinds(j))=0;
      end

      ttt=workdb(curconf).prop.sdesc(end);
      if ttt=='S' || ttt=='T' || ttt=='V'
         conf_syn(i)=conf_syn(i)+1;
      elseif ttt=='A' || ttt=='B' || ttt=='C'
         conf_anti(i)=conf_anti(i)+1;
      else
         error(['Unknown conformation chintype detected: ' workdb(curconf).prop.sdesc]);
      end

      if workdb(curconf).prop.Pdeg<90 || workdb(curconf).prop.Pdeg>270
         conf_north(i)=conf_north(i)+1;
      else
         conf_south(i)=conf_south(i)+1;
      end
      if ceil(workdb(curconf).prop.tbeta/120)==1
         conf_betagp(i)=conf_betagp(i)+1;
      elseif ceil(workdb(curconf).prop.tbeta/120)==0
         conf_betagm(i)=conf_betagm(i)+1;
      else
         conf_betat(i)=conf_betat(i)+1;
      end
      if ceil(workdb(curconf).prop.tgamma/120)==1
         conf_gammagp(i)=conf_gammagp(i)+1;
      elseif ceil(workdb(curconf).prop.tgamma/120)==0
         conf_gammagm(i)=conf_gammagm(i)+1;
      else
         conf_gammat(i)=conf_gammat(i)+1;
      end
      if ceil(workdb(curconf).prop.tdelta/120)==1
         conf_deltagp(i)=conf_deltagp(i)+1;
      elseif ceil(workdb(curconf).prop.tdelta/120)==0
         conf_deltagm(i)=conf_deltagm(i)+1;
      else
         conf_deltat(i)=conf_deltat(i)+1;
      end
      if ceil(workdb(curconf).prop.tepsilon/120)==1
         conf_epsilongp(i)=conf_epsilongp(i)+1;
      elseif ceil(workdb(curconf).prop.tepsilon/120)==0
         conf_epsilongm(i)=conf_epsilongm(i)+1;
      else
         conf_epsilont(i)=conf_epsilont(i)+1;
      end
      if (moltype>100 && mod(moltype,100)==40) || moltype==7
          if ceil(workdb(curconf).prop.teta/120)==1
             conf_etagp(i)=conf_etagp(i)+1;
          elseif ceil(workdb(curconf).prop.teta/120)==0
             conf_etagm(i)=conf_etagm(i)+1;
          else
             conf_etat(i)=conf_etat(i)+1;
          end
          if ceil(workdb(curconf).prop.tteta/120)==1
             conf_tetagp(i)=conf_tetagp(i)+1;
          elseif ceil(workdb(curconf).prop.tteta/120)==0
             conf_tetagm(i)=conf_tetagm(i)+1;
          else
             conf_tetat(i)=conf_tetat(i)+1;
          end
      end

      conf_Htype_popul(i)= conf_Htype_popul(i)+population(curconf); 
    
    end
  %if i==7 , keyboard, end
      
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

    if 0
      hold on
      plot(AIM.AHBang(curinds),AIM.HBdist(curinds),[plotsigns(mod(iplotsigns,numel(plotsigns))+1) 'k'])
      xlabel('AH...B angle,degrees');
      ylabel('H...B distance,A');
      title(Hbondtype.title{i});
      iplotsigns=iplotsigns+1;
      keyboard
    end

end

if any(AIM.pinds(find(AIM.ABdist==0),3)), error('Zero ABdist values detected.'), end %analyze only bonds with three nonzero pinds
if any(AIM.HBdist==0), error('Zero HBdist values detected.'), end
if any(AIM.pinds(find(AIM.ABdist==0),3)), error('Zero AHBang values detected.'), end


%min(pC1pH12bO2notAIManti.HBdist)
%max(pC1pH12bO2notAIManti.HBdist)
%min(pC1pH12bO2notAIManti.AHBang)
%max(pC1pH12bO2notAIManti.AHBang)

%-------------------------------------------------------
    if flanalyzeC6H6
        %determining C6H6 bond parameters sensetivity (bond length, VM frequency, VM intensity) to H-bond formation
        mean_bC6bH6pO4pO5notAIManti_AHdist = mean(bC6bH6pO4pO5notAIManti.AHdist) %#ok
        nstd_bC6bH6pO4pO5notAIManti_AHdist = std(bC6bH6pO4pO5notAIManti.AHdist)/mean_bC6bH6pO4pO5notAIManti_AHdist %#ok
        mean_bC6bH6pO4_AIManti_AHdist = mean(bC6bH6pO4_AIManti.AHdist) %#ok
        nstd_bC6bH6pO4_AIManti_AHdist = std(bC6bH6pO4_AIManti.AHdist)/mean_bC6bH6pO4_AIManti_AHdist %#ok
        mean_bC6bH6pO5_AIManti_AHdist = mean(bC6bH6pO5_AIManti.AHdist) %#ok
        nstd_bC6bH6pO5_AIManti_AHdist = std(bC6bH6pO5_AIManti.AHdist)/mean_bC6bH6pO5_AIManti_AHdist %#ok
        mean_bC6bH6pO4pO5_AIManti_AHdist = mean(bC6bH6pO4pO5_AIManti.AHdist) %#ok
        nstd_bC6bH6pO4pO5_AIManti_AHdist = std(bC6bH6pO4pO5_AIManti.AHdist)/mean_bC6bH6pO4pO5_AIManti_AHdist %#ok

        mean_bC6bH6pO4pO5notAIManti_VMfreq = mean(bC6bH6pO4pO5notAIManti.VMfreq) %#ok
        nstd_bC6bH6pO4pO5notAIManti_VMfreq = std(bC6bH6pO4pO5notAIManti.VMfreq)/mean_bC6bH6pO4pO5notAIManti_VMfreq %#ok
        mean_bC6bH6pO4_AIManti_VMfreq = mean(bC6bH6pO4_AIManti.VMfreq) %#ok
        nstd_bC6bH6pO4_AIManti_VMfreq = std(bC6bH6pO4_AIManti.VMfreq)/mean_bC6bH6pO4_AIManti_VMfreq %#ok
        mean_bC6bH6pO5_AIManti_VMfreq = mean(bC6bH6pO5_AIManti.VMfreq) %#ok
        nstd_bC6bH6pO5_AIManti_VMfreq = std(bC6bH6pO5_AIManti.VMfreq)/mean_bC6bH6pO5_AIManti_VMfreq %#ok
        mean_bC6bH6pO4pO5_AIManti_VMfreq = mean(bC6bH6pO4pO5_AIManti.VMfreq) %#ok
        nstd_bC6bH6pO4pO5_AIManti_VMfreq = std(bC6bH6pO4pO5_AIManti.VMfreq)/mean_bC6bH6pO4pO5_AIManti_VMfreq %#ok

        mean_bC6bH6pO4pO5notAIManti_VMint = mean(bC6bH6pO4pO5notAIManti.VMint) %#ok
        nstd_bC6bH6pO4pO5notAIManti_VMint = std(bC6bH6pO4pO5notAIManti.VMint)/mean_bC6bH6pO4pO5notAIManti_VMint %#ok
        mean_bC6bH6pO4_AIManti_VMint = mean(bC6bH6pO4_AIManti.VMint) %#ok
        nstd_bC6bH6pO4_AIManti_VMint = std(bC6bH6pO4_AIManti.VMint)/mean_bC6bH6pO4_AIManti_VMint %#ok
        mean_bC6bH6pO5_AIManti_VMint = mean(bC6bH6pO5_AIManti.VMint) %#ok
        nstd_bC6bH6pO5_AIManti_VMint = std(bC6bH6pO5_AIManti.VMint)/mean_bC6bH6pO5_AIManti_VMint %#ok
        mean_bC6bH6pO4pO5_AIManti_VMint = mean(bC6bH6pO4pO5_AIManti.VMint) %#ok
        nstd_bC6bH6pO4pO5_AIManti_VMint = std(bC6bH6pO4pO5_AIManti.VMint)/mean_bC6bH6pO4pO5_AIManti_VMint %#ok
        %min_bC6bH6pO4notAIManti_AHBang = min(bC6bH6pO4notAIManti.AHBang)
        %max_bC6bH6pO4notAIManti_AHBang = max(bC6bH6pO4notAIManti.AHBang)
    end

    if flanalyze_pO2HbO2
        %determining C2'H bond parameters sensetivity (bond length, VM frequency, VM intensity) to H-bond formation
        mean_pO2pH22bO2notAIManti_AHdist = mean(pO2pH22bO2notAIManti.AHdist) %#ok
        nstd_pO2pH22bO2notAIManti_AHdist = std(pO2pH22bO2notAIManti.AHdist)/mean_pO2pH22bO2notAIManti_AHdist %#ok
        mean_pO2pH22bO2_AIManti_AHdist = mean(pO2pH22bO2_AIManti.AHdist) %#ok
        nstd_pO2pH22bO2_AIManti_AHdist = std(pO2pH22bO2_AIManti.AHdist)/mean_pO2pH22bO2_AIManti_AHdist %#ok

        mean_pO2pH22bO2notAIManti_VMfreq = mean(pO2pH22bO2notAIManti.VMfreq) %#ok
        nstd_pO2pH22bO2notAIManti_VMfreq = std(pO2pH22bO2notAIManti.VMfreq)/mean_pO2pH22bO2notAIManti_VMfreq %#ok
        mean_pO2pH22bO2_AIManti_VMfreq = mean(pO2pH22bO2_AIManti.VMfreq) %#ok
        nstd_pO2pH22bO2_AIManti_VMfreq = std(pO2pH22bO2_AIManti.VMfreq)/mean_pO2pH22bO2_AIManti_VMfreq %#ok

        mean_pO2pH22bO2notAIManti_VMint = mean(pO2pH22bO2notAIManti.VMint) %#ok
        nstd_pO2pH22bO2notAIManti_VMint = std(pO2pH22bO2notAIManti.VMint)/mean_pO2pH22bO2notAIManti_VMint %#ok
        mean_pO2pH22bO2_AIManti_VMint = mean(pO2pH22bO2_AIManti.VMint) %#ok
        nstd_pO2pH22bO2_AIManti_VMint = std(pO2pH22bO2_AIManti.VMint)/mean_pO2pH22bO2_AIManti_VMint %#ok
    end
    if flanalyze_pO5HbO2
        %determining C2'H bond parameters sensetivity (bond length, VM frequency, VM intensity) to H-bond formation
        mean_pO5pH53bO2notAIMsyn_AHdist = mean(pO5pH53bO2notAIMsyn.AHdist) %#ok
        nstd_pO5pH53bO2notAIMsyn_AHdist = std(pO5pH53bO2notAIMsyn.AHdist)/mean_pO5pH53bO2notAIMsyn_AHdist %#ok
        mean_pO5pH53bO2_AIMsyn_AHdist = mean(pO5pH53bO2_AIMsyn.AHdist) %#ok
        nstd_pO5pH53bO2_AIMsyn_AHdist = std(pO5pH53bO2_AIMsyn.AHdist)/mean_pO5pH53bO2_AIMsyn_AHdist %#ok

        mean_pO5pH53bO2notAIMsyn_VMfreq = mean(pO5pH53bO2notAIMsyn.VMfreq) %#ok
        nstd_pO5pH53bO2notAIMsyn_VMfreq = std(pO5pH53bO2notAIMsyn.VMfreq)/mean_pO5pH53bO2notAIMsyn_VMfreq %#ok
        mean_pO5pH53bO2_AIMsyn_VMfreq = mean(pO5pH53bO2_AIMsyn.VMfreq) %#ok
        nstd_pO5pH53bO2_AIMsyn_VMfreq = std(pO5pH53bO2_AIMsyn.VMfreq)/mean_pO5pH53bO2_AIMsyn_VMfreq %#ok

        mean_pO5pH53bO2notAIMsyn_VMint = mean(pO5pH53bO2notAIMsyn.VMint) %#ok
        nstd_pO5pH53bO2notAIMsyn_VMint = std(pO5pH53bO2notAIMsyn.VMint)/mean_pO5pH53bO2notAIMsyn_VMint %#ok
        mean_pO5pH53bO2_AIMsyn_VMint = mean(pO5pH53bO2_AIMsyn.VMint) %#ok
        nstd_pO5pH53bO2_AIMsyn_VMint = std(pO5pH53bO2_AIMsyn.VMint)/mean_pO5pH53bO2_AIMsyn_VMint %#ok
    end
    if flanalyze_pO2HbN3
        %determining C2'H bond parameters sensetivity (bond length, VM frequency, VM intensity) to H-bond formation
        mean_pO2pH22bN3notAIManti_AHdist = mean(pO2pH22bN3notAIManti.AHdist) %#ok
        nstd_pO2pH22bN3notAIManti_AHdist = std(pO2pH22bN3notAIManti.AHdist)/mean_pO2pH22bN3notAIManti_AHdist %#ok
        mean_pO2pH22bN3_AIManti_AHdist = mean(pO2pH22bN3_AIManti.AHdist) %#ok
        nstd_pO2pH22bN3_AIManti_AHdist = std(pO2pH22bN3_AIManti.AHdist)/mean_pO2pH22bN3_AIManti_AHdist %#ok

        pO2pH22bN3notAIManti.sdesc
        pO2pH22bN3_AIManti.sdesc
        pO2pH22bN3notAIManti.VMfreq
        pO2pH22bN3_AIManti.VMfreq
        mean_pO2pH22bN3notAIManti_VMfreq = mean(pO2pH22bN3notAIManti.VMfreq) %#ok
        nstd_pO2pH22bN3notAIManti_VMfreq = std(pO2pH22bN3notAIManti.VMfreq)/mean_pO2pH22bN3notAIManti_VMfreq %#ok
        mean_pO2pH22bN3_AIManti_VMfreq = mean(pO2pH22bN3_AIManti.VMfreq) %#ok
        nstd_pO2pH22bN3_AIManti_VMfreq = std(pO2pH22bN3_AIManti.VMfreq)/mean_pO2pH22bN3_AIManti_VMfreq %#ok

        pO2pH22bN3notAIManti.VMint
        pO2pH22bN3_AIManti.VMint
        mean_pO2pH22bN3notAIManti_VMint = mean(pO2pH22bN3notAIManti.VMint) %#ok
        nstd_pO2pH22bN3notAIManti_VMint = std(pO2pH22bN3notAIManti.VMint)/mean_pO2pH22bN3notAIManti_VMint %#ok
        mean_pO2pH22bN3_AIManti_VMint = mean(pO2pH22bN3_AIManti.VMint) %#ok
        nstd_pO2pH22bN3_AIManti_VMint = std(pO2pH22bN3_AIManti.VMint)/mean_pO2pH22bN3_AIManti_VMint %#ok
    end
    if flanalyze_pO5HbN3
        %determining C2'H bond parameters sensetivity (bond length, VM frequency, VM intensity) to H-bond formation
        mean_pO5pH53bN3notAIMsyn_AHdist = mean(pO5pH53bN3notAIMsyn.AHdist) %#ok
        nstd_pO5pH53bN3notAIMsyn_AHdist = std(pO5pH53bN3notAIMsyn.AHdist)/mean_pO5pH53bN3notAIMsyn_AHdist %#ok
        mean_pO5pH53bN3_AIMsyn_AHdist = mean(pO5pH53bN3_AIMsyn.AHdist) %#ok
        nstd_pO5pH53bN3_AIMsyn_AHdist = std(pO5pH53bN3_AIMsyn.AHdist)/mean_pO5pH53bN3_AIMsyn_AHdist %#ok

        pO5pH53bN3notAIMsyn.sdesc
        pO5pH53bN3_AIMsyn.sdesc
        pO5pH53bN3notAIMsyn.VMfreq
        pO5pH53bN3_AIMsyn.VMfreq
        mean_pO5pH53bN3notAIMsyn_VMfreq = mean(pO5pH53bN3notAIMsyn.VMfreq) %#ok
        nstd_pO5pH53bN3notAIMsyn_VMfreq = std(pO5pH53bN3notAIMsyn.VMfreq)/mean_pO5pH53bN3notAIMsyn_VMfreq %#ok
        mean_pO5pH53bN3_AIMsyn_VMfreq = mean(pO5pH53bN3_AIMsyn.VMfreq) %#ok
        nstd_pO5pH53bN3_AIMsyn_VMfreq = std(pO5pH53bN3_AIMsyn.VMfreq)/mean_pO5pH53bN3_AIMsyn_VMfreq %#ok

        pO5pH53bN3notAIMsyn.VMint
        pO5pH53bN3_AIMsyn.VMint
        mean_pO5pH53bN3notAIMsyn_VMint = mean(pO5pH53bN3notAIMsyn.VMint) %#ok
        nstd_pO5pH53bN3notAIMsyn_VMint = std(pO5pH53bN3notAIMsyn.VMint)/mean_pO5pH53bN3notAIMsyn_VMint %#ok
        mean_pO5pH53bN3_AIMsyn_VMint = mean(pO5pH53bN3_AIMsyn.VMint) %#ok
        nstd_pO5pH53bN3_AIMsyn_VMint = std(pO5pH53bN3_AIMsyn.VMint)/mean_pO5pH53bN3_AIMsyn_VMint %#ok
    end
    
    if flanalyzeC8H8
        %determining C8H8 bond parameters sensetivity (bond length, VM frequency, VM intensity) to H-bond formation
        mean_bC8bH8pO5notAIManti_AHdist = mean(bC8bH8pO5notAIManti.AHdist) %#ok
        nstd_bC8bH8pO5notAIManti_AHdist = std(bC8bH8pO5notAIManti.AHdist)/mean_bC8bH8pO5notAIManti_AHdist %#ok
        mean_bC8bH8pO5_AIManti_AHdist = mean(bC8bH8pO5_AIManti.AHdist) %#ok
        nstd_bC8bH8pO5_AIManti_AHdist = std(bC8bH8pO5_AIManti.AHdist)/mean_bC8bH8pO5_AIManti_AHdist %#ok

        bC8bH8pO5notAIManti.sdesc
        bC8bH8pO5_AIManti.sdesc
        bC8bH8pO5notAIManti.VMfreq
        bC8bH8pO5_AIManti.VMfreq
        mean_bC8bH8pO5notAIManti_VMfreq = mean(bC8bH8pO5notAIManti.VMfreq) %#ok
        nstd_bC8bH8pO5notAIManti_VMfreq = std(bC8bH8pO5notAIManti.VMfreq)/mean_bC8bH8pO5notAIManti_VMfreq %#ok
        mean_bC8bH8pO5_AIManti_VMfreq = mean(bC8bH8pO5_AIManti.VMfreq) %#ok
        nstd_bC8bH8pO5_AIManti_VMfreq = std(bC8bH8pO5_AIManti.VMfreq)/mean_bC8bH8pO5_AIManti_VMfreq %#ok

        bC8bH8pO5notAIManti.VMint
        bC8bH8pO5_AIManti.VMint
        mean_bC8bH8pO5notAIManti_VMint = mean(bC8bH8pO5notAIManti.VMint) %#ok
        nstd_bC8bH8pO5notAIManti_VMint = std(bC8bH8pO5notAIManti.VMint)/mean_bC8bH8pO5notAIManti_VMint %#ok
        mean_bC8bH8pO5_AIManti_VMint = mean(bC8bH8pO5_AIManti.VMint) %#ok
        nstd_bC8bH8pO5_AIManti_VMint = std(bC8bH8pO5_AIManti.VMint)/mean_bC8bH8pO5_AIManti_VMint %#ok
    end
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


    xlswrite1spec(File,'Data for all Hbonds detected by AIM','AIMavgdata:',['A' int2str(l+1)],params);
    xlswrite1spec(File,corrdatadesc,'AIMavgdata',['A' int2str(l+2)],params);
    xlswrite1spec(File,corrdata','AIMavgdata',['B' int2str(l+2)],params);

    invoke(Excel.ActiveWorkbook,'Save');
    Excel.Quit
    Excel.delete
    clear Excel
    
  catch
    Excel.Quit
    Excel.delete
    clear Excel
  end

end

toc
diary off
