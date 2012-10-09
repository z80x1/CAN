%%plot energy distribution (energy ladder)
%
% Version 1.0    
% Last modified  R O Zhurakivsky 2008-05-19
% Created        R O Zhurakivsky 2006-?-?

format compact
clear
pindsdef
atomsind

if 1

  if 0
    workdbname=[CD.dbdir filesep 'r412_g_dftV4_or.mat']
    T0=NaN;
    T1=298.15;
    T2=NaN;
  elseif 1
    workdbname=[CD.dbdir filesep 'r240_g_dftV2_or.mat']
    %workdbname=[CD.dbdir filesep 'r12_g_dft420_or.mat'] %#ok
    T0=0;
    T1=298.15;
    T2=420;
  end

  load(workdbname,'workdb')

  recnum=numel(workdb);
  if ~recnum
    error('Database is empty')
  end

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

  desc='';
  energy0=[]; %electron energy
  energyT0=[];%electron energy + ZPE
  energyT1=[];%electron energy + GEC(T1)
  energyT2=[];%electron energy + GEC(T2)

  for i=1:recnum
    if workdb(i).new=='Y'

      ms0=workdb(i);

      desc(end+1,:) = ms0.prop.sdesc;


      if isfield(ms0.gaussian,'T')
        if ~isnan(T1)
            tind1 = find(ms0.gaussian.T==T1);
            if isempty(tind1), error(['Desired temperature T=' num2str(T1) 'K is not found']), end  
        end
        if ~isnan(T2)
            tind2 = find(ms0.gaussian.T==T2);
            if isempty(tind2), error(['Desired temperature T=' num2str(T2) 'K is not found']), end  
        end
      else %???
        if ~isnan(T1)
            if T1~=298.15, warning('Check if correct temperature specified and correct datafile used'), end; 
            tind1=1;
        end
        if ~isnan(T2)
            if T2~=298.15, warning('Check if correct temperature specified and correct datafile used'), end; 
            tind2=1;
        end
      end

      if isfield(ms0.gaussian,energyfield)
        energy0(end+1) = ms0.gaussian.(energyfield);
      else
        energy0(end+1)=inf;
        disp(['warning: #' int2str(i) ': energy field ' energyfield ' not found.']);
      end


      if ~isnan(T0)
         energyT0(end+1) = energy0(end)+ms0.gaussian.ZPE/CC.encoef;
      end
      if ~isnan(T1)
        energyT1(end+1) = energy0(end)+ms0.gaussian.GEC(tind1)/CC.encoef;
      end
      if ~isnan(T2)
          energyT2(end+1) = energy0(end)+ms0.gaussian.GEC(tind2)/CC.encoef;
      end

    end
  end


  [energyT0,ind]=sort(energyT0);
  desc=desc(ind,:);

  energyT0=(energyT0-min(energyT0))*CC.encoef;

  x0=zeros(size(energyT0));
  x1=ones(size(energyT0));

  subplot(1,3,1); %0K
  set(gca,'FontSize',14);
  plot([x0;x1],[energyT0;energyT0],'k')
  title('T=0K','FontSize',16)
  ylabel('\Delta G,kcal/mol','FontSize',16);
  axissaved=axis;
  axissaved=[axissaved(1) axissaved(2) 0 1.05*max(energyT0)];
  axis(axissaved);

  %[energy,ind]=sort(energy);
  %desc=desc(ind,:);

  energyT1=(energyT1-min(energyT1))*CC.encoef;

  subplot(1,3,2); %T1
  set(gca,'FontSize',14);
  en = sort(energyT1);
  plot([x0;x1],[en;en],'k')
  title(['T=' num2str(T1) 'K'],'FontSize',16)
  axis(axissaved);


  energyT2=(energyT2-min(energyT2))*CC.encoef;

  subplot(1,3,3); %T2
  set(gca,'FontSize',14);
  en = sort(energyT2);
  plot([x0;x1],[en;en],'k')
  title(['T=' num2str(T2) 'K'],'FontSize',16)
  axis(axissaved);

end

