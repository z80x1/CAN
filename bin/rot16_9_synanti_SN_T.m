%plot temperature dependence of relavive populations of Syn/Anti and South/North conformers
%and population
%
% Version 1.0    
% Last modified  R O Zhurakivsky 2008-08-01
% Created        R O Zhurakivsky 2006-?-?

clear

atomsind

format compact

  moltype=161  %#ok
  theory='dftV3'  %#ok
%  theory='mp2V2'  %#ok
  usedpackage='Gaussian'  %#ok
  onlyoriginal=1;  % process db with only original conformations
  dbsuffix='' %#ok
%  T=[298.15 320 340 360 380 400 420] %#ok
  T=[298.15 300 320 340 360 380 400 420] %#ok

color='r';
color2='b';

fl_synanti = 1
fl_SN = 1

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
recnum=numel(workdb);
if ~recnum
  error('Database is empty')
end

energy=zeros(recnum,numel(T));
desc=zeros(recnum,numel(workdb(1).prop.sdesc));
[tchi,Pdeg]=deal(zeros(recnum,1));

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
population=repmat(normcoef,recnum,1).*exp(-energy/CC.k./repmat(T,recnum,1));

cla

ylabelChistr='';
if fl_synanti
    indA = find(desc(1:i,end)=='A');
    indS = find(desc(1:i,end)=='S');
    popA = sum(population(indA,:),1);
    popS = sum(population(indS,:),1);

    hold on

    ratioChi=popA./popS;
    ylabelChistr= 'Anti:Syn';
    titleChistr = [ylabelChistr ' equlibrium temperature dependence'];
    
    plot(T,ratioChi,['k' 's-'],'MarkerEdgeColor','Black','MarkerFaceColor',color,'MarkerSize',10,'LineWidth',2);

%    grid on
%    set(gca,'XMinorGrid','on');
%    title(titleChistr,'FontSize',20);
%    xlabel('T, K','FontSize',16);
%    ylabel(ylabelChistr,'FontSize',16);

    axis([ 0.98*min(T) 1.02*max(T) 0.98*min(ratioChi) 1.02*max(ratioChi)]);
    set(gca,'FontSize',16)
end
if fl_SN
    sugarconf=desc(:,1);
%    indC3endo = find(sugarconf=='A');
%    indC4exo =  find(sugarconf=='B');
%    indO4endo = find(sugarconf=='C');
%    indC1exo =  find(sugarconf=='D');
%    indC2endo=  find(sugarconf=='E');
%    indC3exo =  find(sugarconf=='F');
%    indC4endo = find(sugarconf=='G');
%    indC2exo =  find(sugarconf=='J');
%    indNorth = [indC3endo; indC4exo; indC2exo];
%    indSouth = [indC1exo; indC2endo; indC3exo; indC4endo];
    indNorth = find(Pdeg<=90 | Pdeg>=270);
    indSouth = find(abs(Pdeg-180)<90);

    if numel(indNorth)+numel(indSouth)~=recnum
        error('All conformers are not separated by sugar conformation');
    end
    popNorth = sum(population(indNorth,:),1);
    popSouth = sum(population(indSouth,:),1);

    hold on

%    ratioP=popNorth./popSouth;
%    ylabelPstr= 'North:South';
    ratioP=popNorth./popSouth;
    ylabelPstr= 'North:South';
    titlePstr = [ylabelPstr ' equlibrium temperature dependence'];
   

    plot(T,ratioP,['k' 'd-'],'MarkerEdgeColor','Black','MarkerFaceColor',color2,'MarkerSize',10,'LineWidth',2);

%    grid on
%    set(gca,'XMinorGrid','on');
%    title(titlePstr,'FontSize',20);
%    xlabel('T, K','FontSize',16);
%    ylabel(ylabelPstr,'FontSize',16);

    oldaxis=axis;
    newaxis=[ 0.98*min(T) 1.02*max(T) 0.98*min(ratioP) 1.02*max(ratioP)];

    axis([min(oldaxis(1),newaxis(1)) max(oldaxis(2),newaxis(2)) min(oldaxis(3),newaxis(3)) max(oldaxis(4),newaxis(4))]);
%    axis([newaxis(1) newaxis(2) 0.1 0.5]); %rUrd
%    axis([newaxis(1) newaxis(2) 0.5 3.5]); %dUrd dThd
%    axis([newaxis(1) newaxis(2) 0.0 0.35]); %rGuo

    set(gca,'FontSize',16)

    grid on
    set(gca,'XMinorGrid','on');
    xlabel('T, K','FontSize',16);
    ylabel('populations ratio','FontSize',16);
    title(['r' int2str(moltype) ': Populations ratio temperature dependence'],'FontSize',20);
    if isempty(ylabelChistr)
        obj=legend(ylabelPstr,'Location','Best');
    elseif isempty(ylabelPstr)
        obj=legend(ylabelChistr,'Location','Best');
    else
        obj=legend(ylabelChistr,ylabelPstr,'Location','Best');
    end
    set(gca,'XMinorGrid','off')
    set(obj,'FontSize',16)
end

disp('Done!')
