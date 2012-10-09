%%plot energy distribution (energy ladder) for severat molecules and one temperature
%
% Version 1.0    
% Last modified  R O Zhurakivsky 2012-05-04
% Created        R O Zhurakivsky 2012-05-04

format compact
clear
pindsdef
atomsind


  T=298.15
  workdbnames=[ {[CD.dbdir filesep 'r412_g_dftV4_or.mat']};
              ]
  energyfields=[{'MP2_6311__Gdp'};
                ];
  titles=[{'d4T'}];
  DNAlike=[{'CabC'}];
%   workdbnames=[ {[CD.dbdir filesep 'r412_g_dftV4_or.mat']};
%                 {[CD.dbdir filesep 'r512_g_dftV4_or.mat']};
%                 {[CD.dbdir filesep 'r212_g_dftV4_or.mat']};
%                 {[CD.dbdir filesep 'r112_g_dftV3bis_or.mat']};
%                 {[CD.dbdir filesep 'r312_g_dftV3bis_or.mat']};
%               ]
%   energyfields=[{'MP2_6311__Gdp'};
%                 {'MP2_6311__Gdp'};
%                 {'MP2_6311__Gdp'};
%                 {'B3LYP_631Gdp'}; %there are no SP !!!
%                 {'B3LYP_631Gdp'}; %there are no SP !!!
%                 ];
%   titles=[{'d4T'};{'d4U'};{'d4C'};{'d4A'};{'d4G'}];
%   DNAlike=[{'CabC'};{'CabC'};{'CabC'};{'BabC'};{'BabC'}];

  for imol=1:numel(workdbnames)

      load(workdbnames{imol},'workdb')

      recnum=numel(workdb);
      if ~recnum
        error('Database is empty')
      end

      energyfield = energyfields{imol};
      disp(['energy field ' energyfield ' used']);

      desc='';
      energy0=[]; %electron energy
      energyT=[]; 

      iconf = 0;
      for i=1:recnum
        if workdb(i).new=='Y'
          iconf = iconf+1;

          ms0=workdb(i);

          desc(iconf,:) = ms0.prop.sdesc;


          if isfield(ms0.gaussian,'T')
            if ~isnan(T)
                tind = find(ms0.gaussian.T==T);
                if isempty(tind), error(['Desired temperature T=' num2str(T) 'K is not found']), end  
            end
          else %???
            if ~isnan(T)
                if T~=298.15, warning('Check if correct temperature specified and correct datafile used'), end; 
                tind=1;
            end
          end

          if isfield(ms0.gaussian,energyfield)
            energy0(iconf) = ms0.gaussian.(energyfield);
          else
            energy0(iconf)=inf;
            disp(['warning: #' int2str(i) ': energy field ' energyfield ' not found.']);
          end


          energyT(iconf) = energy0(iconf)+ms0.gaussian.GEC(tind)/CC.encoef;

        end
      end



      [energyT,ind]=sort(energyT);
      desc=desc(ind,:);

      energyT=(energyT-min(energyT))*CC.encoef;

      x0=zeros(size(energyT));
      x1=ones(size(energyT));

      subplot(1,numel(workdbnames),imol); 
      hold on
      set(gca,'FontSize',14);
      plot([x0;x1],[energyT;energyT],'k','LineWidth',1.5);

      DNAind = strcmpar(desc,DNAlike{imol});
      if numel(DNAind)>0
        plot([zeros(size(DNAind));ones(size(DNAind))],[energyT(DNAind);energyT(DNAind)],'r','LineWidth',1.5);
        text(1,energyT(DNAind),' \leftarrow','FontSize',16)
      end

      title(titles{imol},'FontSize',16)
      if imol==1
          ylabel('\Delta G,kcal/mol','FontSize',16);
      end
      set(gca,'XTick',[],'Box','on');
      axissaved=axis;
      axissaved=[axissaved(1) axissaved(2) 0 1.05*max(energyT)];
      axis(axissaved);

end
%  text(0.9,0.5,['T=' num2str(T) 'K'],'FontSize',16)

