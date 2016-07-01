%calculate specific property for all conformations (rglyc_avg)
%and population distributions (syn, anti, S, N, C2endo, C3endo)
%
% Version 1.0    
% Last modified  R O Zhurakivsky 2008-05-19
% Created        R O Zhurakivsky 2006-?-?

clear

atomsind

format compact

  moltype=12  %#ok
  theory='dftV2'  %#ok
  usedpackage='Gaussian'  %#ok
  onlyoriginal=1;  % process db with only original conformations
  T1=298.15 %#ok
  T2=420 %#ok

  color='b' %#ok
  color2='r' %#ok

  includelist.do=0    %#ok
  includelist.type='file';
  includelist.name=['includelist'];
  includelist  %#ok

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
workdbname=[CD.dbdir filesep workdbname '.mat'] %#ok

load(workdbname,'workdb')
recnum=numel(workdb);
if ~recnum
  error('Database is empty')
end

rglyc=zeros(0);
rC2H2=zeros(0);
desc='';
energyT1=[];
energyT2=[];
rCNamino=[];
tau03=[];
basecharge=[];
[tau05,tau06]=deal([]);


inclset={};
if includelist.do
if ~isempty(includelist.name)

  if strcmp(includelist.type,'file')

    [fid,emessage]=fopen(includelist.name);
    if fid==-1
      error(['No file ' includelist.name ' with descriptions to include found: ' emessage]);
    end

    tline=fgetl(fid);
	i=1;
    while tline~=-1
  	  inclfiles{i}=tline;
      tline=fgetl(fid);
	  i=i+1;
    end	
    fclose(fid);
    numfiles=i-1;

  end

  if ~numfiles
	disp('No descriptions found to include.')
  else
  end
end
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

baseatoms=[];
for i=1:numel(pind.labels)
  if pind.labels{i}(1)=='b', baseatoms(end+1)=i; end;
end

for i=1:recnum
    if workdb(i).new=='Y'

        ms0=workdb(i);

        iC1 = ms0.ind(find(find(strcmp(pind.labels,'pC1'))==ms0.pind)); 
        ibN1= ms0.ind(find(find(strcmp(pind.labels,'bN1'))==ms0.pind)); 
        iC2 = ms0.ind(find(find(strcmp(pind.labels,'pC2'))==ms0.pind)); 
        iH22= ms0.ind(find(find(strcmp(pind.labels,'pH22'))==ms0.pind)); 

        desc(end+1,:) = ms0.prop.sdesc;

        rglyc(end+1)=adist(ms0,iC1,ibN1);
        rC2H2(end+1)=adist(ms0,iC2,iH22);

        tau05(end+1)=ms0.prop.tau05;
        tau06(end+1)=ms0.prop.tau06;

    %----------------------
       if isfield(ms0.gaussian,'T')
         tind1 = find(ms0.gaussian.T==T1);
         if isempty(tind1), error(['Desired temperature T=' num2str(T1) 'K is not found']), end  
         tind2 = find(ms0.gaussian.T==T2);
         if isempty(tind2), error(['Desired temperature T=' num2str(T2) 'K is not found']), end  
       else %???
         if T1~=298.15, warning('Check if correct temperature specified and correct datafile used'), end; 
         if T2~=298.15, warning('Check if correct temperature specified and correct datafile used'), end; 
         tind1=1;
         tind2=1;
       end

       if isfield(ms0.gaussian,energyfield)
      	energyT1(end+1) = ms0.gaussian.(energyfield)  + ms0.gaussian.GEC(tind1)/CC.encoef;
      	energyT2(end+1) = ms0.gaussian.(energyfield)  + ms0.gaussian.GEC(tind2)/CC.encoef;
       else
        energyT1(end+1)=inf;
        energyT2(end+1)=inf;
    	disp(['warning: #' int2str(i) ': energy field ' energyfield ' not found.']);
       end


    end

%   [X,iA,iB]=intersect(ms0.pind,baseatoms);
%   basecharge(end+1)=sum(ms0.gaussian.mcharge(iA));


end

rglyc_avg=mean(rglyc) %#ok


%determining conformational equilibrium


%[energy,sortind]=sort(energy);
%desc=desc(sortind,:);
%%!! warning energy and desc are ordered by energy now. !!!

energyT1=(energyT1-min(energyT1))*CC.hartree;
normcoef=1/(sum(exp(-energyT1/CC.k/T1)));
populationT1=normcoef*exp(-energyT1/CC.k/T1);

%[A,IA,IB] = intersect(desc,inclfiles);
%populationT1=normcoef*exp(-energyT1(IA)/CC.k/T1);
%sum(populationT1)
%return

energyT2=(energyT2-min(energyT2))*CC.hartree;
normcoef=1/(sum(exp(-energyT2/CC.k/T2)));
populationT2=normcoef*exp(-energyT2/CC.k/T2);

[population_orderedT1,II] = sort(populationT1,2,'descend');
population_orderedT1 %#ok
desc_orderedT1=desc(II,:);

[population_orderedT2,IIT2] = sort(populationT2,2,'descend');
population_orderedT2 %#ok
desc_orderedT2=desc(IIT2,:);

if 0

    popA = sum(population(find(desc(:,5)=='A'))) %#ok
    popS = sum(population(find(desc(:,5)=='S'))) %#ok

    sugarconf=desc(:,1);
    popNorth = sum(population([find(sugarconf=='A'); find(sugarconf=='B'); find(sugarconf=='I'); find(sugarconf=='J')])) %#ok
    popSouth = sum(population([find(sugarconf=='D'); find(sugarconf=='E'); find(sugarconf=='F'); find(sugarconf=='G')])) %#ok

    popC3endo = sum(population(find(sugarconf=='A'))) %#ok
    popC2endo = sum(population(find(sugarconf=='E'))) %#ok


    drCNamino=rCNamino-min(rCNamino) %#ok
    [drCNamino,I]=sort(drCNamino);
    tau03=tau03(I);
    plot(drCNamino,tau03); 
    title('Correlation between drCNamino and C4N out of NH2 plane angle')
    xlabel('delta rCNamino, A');
    ylabel('tau03,degree');

    [r,p]=corrcoef(drCNamino,tau03) %#ok
end

hold on

if 0

    popsumT1=zeros(1,numel(populationT1)+1);
    for i=0:numel(populationT1)
      popsumT1(i+1)=sum(population_orderedT1(1:i));
    end

    popsumT2=zeros(1,numel(populationT2)+1);
    for i=0:numel(populationT2)
      popsumT2(i+1)=sum(population_orderedT2(1:i));
    end

    plot(0:numel(populationT1),popsumT1,[color 'd-'],'MarkerEdgeColor',color,'MarkerFaceColor',color,'MarkerSize',4);
    plot(0:numel(populationT2),popsumT2,[color2 'o-'],'MarkerEdgeColor',color2,'MarkerFaceColor',color2,'MarkerSize',4);
    
    grid on
    set(gca,'XMinorGrid','on');
    title('Dependence of population of conformers over their number','FontSize',20);
    xlabel('number of conformers','FontSize',16);
    ylabel('relative population','FontSize',16);

    axis([0 recnum 0 1.05]);
    set(gca,'FontSize',16)
    obj=legend(['T=' num2str(T1) 'K'],['T=' num2str(T2) 'K'],'Location','Best');
    set(obj,'FontSize',16)
	
end

if 1
    popAsumT1=zeros(1,numel(population_orderedT1)+1);
    popSsumT1=zeros(1,numel(population_orderedT1)+1);
    for i=0:numel(population_orderedT1)
      indA = find(desc_orderedT1(1:i,5)=='A');
      indS = find(desc_orderedT1(1:i,5)=='S');
      popAsumT1(i+1)=sum(population_orderedT1(indA));
      popSsumT1(i+1)=sum(population_orderedT1(indS));
    end
    popS_popA_limitT1 = popSsumT1(end)./popAsumT1(end);

    popAsumT2=zeros(1,numel(population_orderedT2)+1);
    popSsumT2=zeros(1,numel(population_orderedT2)+1);
    for i=0:numel(population_orderedT2)
      indA = find(desc_orderedT2(1:i,5)=='A');
      indS = find(desc_orderedT2(1:i,5)=='S');
      popAsumT2(i+1)=sum(population_orderedT2(indA));
      popSsumT2(i+1)=sum(population_orderedT2(indS));
    end
    popS_popA_limitT2 = popSsumT2(end)./popAsumT2(end);

    
    
    hold on
%    plot(0:numel(population),popAsum,[color '.-']);
%    plot(0:numel(population),popSsum,[color2 '.-']);
    plot(0:numel(populationT1),popSsumT1./popAsumT1,[color 'o-'],'MarkerEdgeColor',color,'MarkerFaceColor',color,'MarkerSize',4);
    plot(0:numel(populationT2),popSsumT2./popAsumT2,[color2 'd-'],'MarkerEdgeColor',color,'MarkerFaceColor',color,'MarkerSize',4);

    grid on
    set(gca,'XMinorGrid','on');
    title('Dependence of population of Syn/Anti conformers over their number');
    xlabel('number of conformers');
    ylabel('population');

    maxy = max([max(popSsumT1./popAsumT1) max(popSsumT2./popAsumT2)]);
    if maxy==Inf, maxy=2.5;, end;
    axis([ 0 recnum 0 1.15*maxy ]);
    legend(['T=' num2str(T1) 'K'],['T=' num2str(T2) 'K'],'Location','Best');

    plot(0:numel(populationT1), repmat(0.90*popS_popA_limitT1,1,numel(populationT1)+1),'k' );
    plot(0:numel(populationT1), repmat(1.10*popS_popA_limitT1,1,numel(populationT1)+1),'k' );
    plot(0:numel(populationT1), repmat(0.95*popS_popA_limitT1,1,numel(populationT1)+1),'b' );
    plot(0:numel(populationT1), repmat(1.05*popS_popA_limitT1,1,numel(populationT1)+1),'b' );
    plot(0:numel(populationT2), repmat(0.90*popS_popA_limitT2,1,numel(populationT2)+1),'k' );
    plot(0:numel(populationT2), repmat(1.10*popS_popA_limitT2,1,numel(populationT2)+1),'k' );
    plot(0:numel(populationT2), repmat(0.95*popS_popA_limitT2,1,numel(populationT2)+1),'b' );
    plot(0:numel(populationT2), repmat(1.05*popS_popA_limitT2,1,numel(populationT2)+1),'b' );
%    axis([ 0 recnum 0 1.05*max([popAsum(end) popSsum(end)]) ]);
%    legend('Anti','Syn','Location','Best');
    title('Syn/Anti');
	

end

disp('Done!')
