%plot temperature dependence of relavive populations of Syn/Anti and South/North conformers
%and population for all molecules
%NOT COMPLETED!!!
%
% Version 1.0    
% Last modified  R O Zhurakivsky 2010-11-03
% Created        R O Zhurakivsky 2010-11-03

clear

atomsind
format compact

moltypes=[9 12 13 15 16 140 240 340 540]  %#ok
%moltypes=[12 13]  %#ok
theory='dftV2'  %#ok
%  theory='mp2V2'  %#ok
usedpackage='Gaussian'  %#ok
onlyoriginal=1;  % process db with only original conformations
dbsuffix='' %#ok
%  T=[298.15 320 340 360 380 400 420] %#ok
T=[298.15 320 340 360 380 400 420] %#ok
%  conf2proc = 89;%dCyd %use if not all structures have needed data
%  conf2proc = 88;%dAdo %use if not all structures have needed data

color='r'; 
color2='b';

fl_synanti = 1 %plot syn/anti dependence
fl_SN = 0      %plot S/N dependence

for n=1:numel(moltypes)
  moltype = moltypes(n);


  workdbname=['r' int2str(moltype)]   %#ok
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
  if ~isempty(dbsuffix)
      workdbname = [workdbname '_' dbsuffix];
  end
  workdbname=[CD.dbdir filesep workdbname '.mat'] %#ok

  load(workdbname,'workdb')

  if exist('conf2proc','var')
  	recnum = conf2proc;
  else
  	recnum = numel(workdb);
  end
  if ~recnum
    error('Database is empty')
  end

  [tchi,Pdeg]=deal(zeros(recnum,1));
  energy=zeros(recnum,numel(T));
  desc=zeros(recnum,numel(workdb(1).prop.sdesc));

  if isfield(workdb(1).gaussian,'MP2_311__G2dfpd')
   	energyfield='MP2_311__G2dfpd';
  elseif isfield(workdb(1).gaussian,'MP2_6311__G2dfpd')
   	energyfield='MP2_6311__G2dfpd';
  elseif isfield(workdb(1).gaussian,'MP2_6311__Gdp')
   	energyfield='MP2_6311__Gdp';
  else
  	error('energy field not found');
  end
  disp(['energy field ' energyfield ' detected']);

  for i=1:recnum
     
      if workdb(i).new=='Y'

         ms0=workdb(i);
         desc(i,:) = ms0.prop.sdesc;

         tind=zeros(size(T));

         for j=1:numel(T)
            if isfield(ms0.gaussian,'T')
              indT=find(ms0.gaussian.T==T(j));
              if isempty(indT), error(['(#' int2str(i) ') Desired temperature T=' num2str(T(j)) 'K is not found']), end  
              tind(j) = indT;
            else
              if T(j)~=298.15, warning('synanti:corrT','Check if correct temperature specified and correct datafile used'), end; 
              tind(j)=1;
            end

            if isfield(ms0.gaussian,energyfield)
                energy(i,j) = ms0.gaussian.(energyfield)  + ms0.gaussian.GEC(tind(j))/CC.encoef;
            else
                energy(i,j)=inf;
         	    disp(['warning: #' int2str(i) ': energy field ' energyfield ' not found.']);
            end

         end

  %       tchi(i) = ms0.prop.tchi;
         Pdeg(i) = ms0.prop.Pdeg;

      end
  end


  %determining conformational equilibrium

  energymin=min(energy,[],1);
  energy=(energy-repmat(energymin,recnum,1))*CC.hartree;
  normcoef=1./(sum(exp(-energy/CC.k./repmat(T,recnum,1)),1));

  population(n)={repmat(normcoef,recnum,1).*exp(-energy/CC.k./repmat(T,recnum,1))};
  Pdeg_all(n) = {Pdeg};
  desc_all(n) = {desc};

end %n=1:moltypes


cla
hold on
%    grid on
%    set(gca,'XMinorGrid','on');
%    title(titleChistr,'FontSize',20);
%    xlabel('T, K','FontSize',16);
%    ylabel(ylabelChistr,'FontSize',16);
axis([ 0.98*min(T) 1.02*max(T) 0 100]);
set(gca,'FontSize',16)

ylabelChistr='';
for n=1:numel(moltypes)
if fl_synanti
    desc = desc_all{n};
    
    indA = find(desc(:,end)=='A');
    indA = [indA; find(desc(:,end)=='B')];
    indA = [indA; find(desc(:,end)=='C')];
    indS = find(desc(:,end)=='S');
    indS = [indS; find(desc(:,end)=='T')];
    indS = [indS; find(desc(:,end)=='V')];
    popA = sum(population{n}(indA,:),1);
    popS = sum(population{n}(indS,:),1);

    %ratioChi=popA./popS;
    %ylabelChistr= 'Anti:Syn';
    %titleChistr = [ylabelChistr ' equlibrium temperature dependence'];
    
%    plot(T,ratioChi,['k' 's-'],'MarkerEdgeColor','Black','MarkerFaceColor',color,'MarkerSize',10,'LineWidth',2);
    plot(T,100*popA,['k' 's-'],'MarkerEdgeColor','Black','MarkerFaceColor',color,'MarkerSize',10,'LineWidth',2);
    plot(T,100*popS,['k' 's-'],'MarkerEdgeColor','Black','MarkerFaceColor',color2,'MarkerSize',10,'LineWidth',2);

end
end

for n=1:numel(moltypes)
if fl_SN
%    sugarconf=desc(:,1);
    Pdeg = Pdeg_all{n}
    indNorth = find(Pdeg<=90 | Pdeg>=270);
    indSouth = find(abs(Pdeg-180)<90);

    if numel(indNorth)+numel(indSouth)~=recnum
        error('All conformers are not separated by sugar conformation');
    end
    popNorth = sum(population{n}(indNorth,:),1);
    popSouth = sum(population{n}(indSouth,:),1);

%    ratioP=popNorth./popSouth;
%    plot(T,ratioP,['k' 'd-'],'MarkerEdgeColor','Black','MarkerFaceColor',color2,'MarkerSize',10,'LineWidth',2);
    plot(T,100*popNorth,['k' 'd-'],'MarkerEdgeColor','Black','MarkerFaceColor',color,'MarkerSize',10,'LineWidth',2);
    plot(T,100*popSouth,['k' 'd-'],'MarkerEdgeColor','Black','MarkerFaceColor',color2,'MarkerSize',10,'LineWidth',2);

end
end

    grid on
    set(gca,'XMinorGrid','on');
    xlabel('T, K','FontSize',16);
    ylabel('population','FontSize',16);
    title(['r' int2str(moltype) ': Populations temperature dependence'],'FontSize',20);
%     if isempty(ylabelChistr)
%         obj=legend(ylabelPstr,'Location','Best');
%     elseif isempty(ylabelPstr)
%         obj=legend(ylabelChistr,'Location','Best');
%     else
%         obj=legend(ylabelChistr,ylabelPstr,'Location','Best');
%     end
    set(gca,'XMinorGrid','off')
%    set(obj,'FontSize',16)

hold off
disp('Done!')
