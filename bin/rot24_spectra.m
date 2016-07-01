%rot24: make molecule vibrational spectra
%
% Version 1.0    
% Last modified  R O Zhurakivsky 2008-05-19
% Created        R O Zhurakivsky 2006-?-?

%060623: [-] for calculating intencities for reper frequencies all
%frequncies were scaled while repers ones weren't

clear
atomsind
format compact
flPlot=0;
labelHval=0;
labelHval2=0;

%--------------------------------------------------------------
moltype=521

if moltype==12
    molname='dUrd'
    workdbname=[CD.dbdir filesep 'r12_g_dftV2_or.mat'] 
elseif moltype==13
    molname='dThd'; 
    workdbname=[CD.dbdir filesep 'r13_g_dftV2_or.mat'] 
elseif moltype==140
    molname='rAdo'; 
    workdbname=[CD.dbdir filesep 'r140_g_dftV2_or.mat'] 
elseif moltype==521
    molname='d5BrU'; 
    workdbname=[CD.dbdir filesep 'r521_g_dftV2_or.mat'] 
elseif moltype==540
    molname='rUrd'; 
    workdbname=[CD.dbdir filesep 'r540_g_dftV2_or.mat'] 
else
    error('Unknown molecule.');
end

%T=420; %rThd evaporation temperature 
%T=421; %rUrd evaporation temperature 
%T=456; % rAdo evaporation temperature is +183C
T=298.15; 

%freqlimits=[3350 3700]
%freqlimits=[0 4000]

sortmode=2  %#ok %without sorting/sdesc/energy/manual
fl_plotgauss=1 %#ok %gaussians envelopes of integral intensity of modes
fl_plotlines=1 %#ok %plot mode intensity lines

%numconf2plot=[1:111];
%numconf2plot=[6 20 94];
%numconf2plot=[6 20 92];
%numconf2plot=[6 20 137];
numconf2plot=[6 20 111];

%conf2exclude = [{'BaacS'},{'BaaaS'}]; %these conformers will not be included in spectra
conf2exclude = [];

labelHval=[{'pH22'}]  %label of H atom of valence vibrational mode. If specified - only these line in spectra will be plotted in 1/4 plot
labelHval2=[{'pH32'}]  %label of H atom of valence vibrational mode. If specified - only these line in spectra will be plotted in 2/4 plot
fl_plotcurconf=0 %plot first deconvolution section
fl_plotdeconvolsect1=1 %plot first deconvolution section
fl_plotdeconvolsect2=1 %plot second deconvolution section
%--------------------------------------------------------------

clf
hold on
color='kbrgmyc';

if moltype==12
    freqfactor=3428/3617 %experimental to our calc ratio (3428 - dUrd NH3 frequency by Ivanov; 3617 - averaged NH3 freq by our resuts )
elseif moltype==13
    freqfactor=3428/3614 %experimental to our calc ratio (3427.2 - dThd NH3 frequency by Ivanov; 3614 - averaged NH3 freq by our resuts )
%elseif moltype==140
%    freqfactor=1 %experimental to our calc ratio (3427.2 - rUrd NH3 frequency by Ivanov; 3616 - averaged NH3 freq over 111 confs by our resuts )
elseif moltype==540
    freqfactor=3427.2/3616 %experimental to our calc ratio (3427.2 - rUrd NH3 frequency by Ivanov; 3616 - averaged NH3 freq over 111 confs by our resuts )
else
    warning('Cannot find modes frequency scaling data for specified molecule. Defaults used.');
    freqfactor=0.9608;
end


fignum1=100;
fignum2=101;

reperfreq=[3048 3065 3427.2 3492 3568 3614 3633];  %rc frequencies to sum intencities over conformers in 2*halfwidth diapason


load(workdbname,'workdb')


recnum=numel(workdb);
for i=1:recnum

    sdesc(i) = {workdb(i).prop.sdesc};
    if workdb(i).new~='Y'
        energy(i)=inf;
        continue
    end

    if isfield(workdb(i).gaussian,'T')
      tind = find(workdb(i).gaussian.T==T);
      if isempty(tind), error(['Desired temperature T=' num2str(T) 'K is not found']), end  
    else %???
      tind=1;
    end

    if isfield(workdb(i).gaussian,'MP2_6311__Gdp') & isfield(workdb(i).gaussian,'GEC')
        energy(i) = workdb(i).gaussian.MP2_6311__Gdp + workdb(i).gaussian.GEC(tind)/CC.encoef;
    else
        disp('needed energy fields are not found');
        sortbyenergy=0;
    end
end
minenergy=min(energy); %min of E_SP+GEC

if sortmode==0 %without sorting
  ind=1:recnum;
elseif sortmode==1 %sort by sdesc
  [xxx,ind]=sort(sdesc);
elseif sortmode==2 %sort by energy
  [xxx,ind]=sort(energy);
elseif sortmode==3 %select
  ss=['AcbaA';'BcbaA';'AccaA';'AcbaS';'BccaS'];
  ind=[];   
  for i=1:size(ss,1)
    ind(end+1)=strcmpcellar(sdesc,ss(i,:))
  end

else 
  error('invalid sort mode');
end
confnumber='';
j=1; %calculating real number of cycles over i

N=4000; %number of points to calculate spectra in
%if freqlimits(1)>freqlimits(2)
%   error('first element of freqlimits must be less than the second');
%end
%N=ceil(freqlimits(2)-freqlimits(1));
sumY=zeros(3,N);

reperint =  zeros(numel(reperfreq),1);
reperintnorm =  zeros(numel(reperfreq),recnum);

freqmodenum = numel(workdb(1).freq.freq);

for ijk=1:numel(labelHval)
  pindHval(ijk)=strcmpcellar(pind.labels,labelHval{ijk});
end
for ijk=1:numel(labelHval2)
  pindHval2(ijk)=strcmpcellar(pind.labels,labelHval2{ijk});
end

ynormsum{1}=zeros(recnum,freqmodenum); %sum over mode intencities is now calculated only for one mode
ynormsum{2}=zeros(recnum,freqmodenum); %sum over mode intencities is now calculated only for one mode
ynormsum{3}=zeros(recnum,freqmodenum); %sum over mode intencities is now calculated only for one mode

hgauss1=[];
hgauss2=[];
hgauss3=[];
hgauss4=[];

for i=ind %sorting conformers
%disp(['Conformer #' int2str(j)]);

    if ~isempty(confnumber)
      i=confnumber; %#ok
    end

    ms0=workdb(i);

    if strcmpcellar(conf2exclude,ms0.prop.sdesc)
        continue
    end

    if isfield(ms0.gaussian,'T')
      tind = find(ms0.gaussian.T==T);
      if isempty(tind), error(['Desired temperature T=' num2str(T) 'K is not found']), end  
    else %???
      tind=1;
    end

    freqscaled = workdb(i).freq.freq*freqfactor;

%    if pindHval
      [XXX,indVMarr,III]=intersect(ms0.pind(strcmpcellar(ms0.labels,'H')),pindHval);
      modeind{2} = ms0.prop.HbondVMind(indVMarr);
      [XXX,indVMarr,III]=intersect(ms0.pind(strcmpcellar(ms0.labels,'H')),pindHval2);
      modeind{3} = ms0.prop.HbondVMind(indVMarr);
%    else     

%      modeind{1} = find(freqscaled>freqlimits(1) & freqscaled<freqlimits(2));
      modeind{1} = 1:freqmodenum;
%    end

    x{1}=freqscaled(modeind{1});
    x{2}=freqscaled(modeind{2});
    x{3}=freqscaled(modeind{3});
    if isfield(ms0.freq,'intencity')
        y{1}=ms0.freq.intencity(modeind{1});
        y{2}=ms0.freq.intencity(modeind{2});
        y{3}=ms0.freq.intencity(modeind{3});
  %     ynorm=y*exp((minenergy-ms0.GO.energy)*CC.hartree/CC.k/T);
        ynorm{1}=y{1}*exp((minenergy-workdb(i).gaussian.MP2_6311__Gdp - workdb(i).gaussian.GEC(tind)/CC.encoef)*CC.hartree/CC.k/T);
        ynorm{2}=y{2}*exp((minenergy-workdb(i).gaussian.MP2_6311__Gdp - workdb(i).gaussian.GEC(tind)/CC.encoef)*CC.hartree/CC.k/T);
        ynorm{3}=y{3}*exp((minenergy-workdb(i).gaussian.MP2_6311__Gdp - workdb(i).gaussian.GEC(tind)/CC.encoef)*CC.hartree/CC.k/T);
%        ynorm=y*exp((minenergy-energy(i)/CC.encoef)*CC.hartree/CC.k/T);
    else
        y=zeros(size(x{1}));
        ynorm=zeros(size(x{1}));
    end
    ynormsum{1}(j+1,1:numel(modeind{1}))=ynormsum{1}(j,1:numel(modeind{1}))+ynorm{1};
    ynormsum{2}(j+1,1:numel(modeind{2}))=ynormsum{2}(j,1:numel(modeind{2}))+ynorm{2};
    ynormsum{3}(j+1,1:numel(modeind{3}))=ynormsum{3}(j,1:numel(modeind{3}))+ynorm{3};

    if j==1
        xlen=1.15*max(x{1});
        X=0:xlen/(N-1):xlen;

%        xlen(1)=0.95*freqlimits(1);
%        xlen(2)=1.05*freqlimits(2);
%        X=xlen(1):(xlen(2)-xlen(1))/(N-1):xlen(2);
    end

    %a=1/5^2;
    %halfwidth=2*sqrt(1/a)*sqrt(log(2))  - halfwidth of gaussian

    halfwidth=8.5; %cm^-1
  %  halfwidth=10; %cm^-1  for dThd (bad)

    a=4*log(2)/halfwidth^2;

    for I=1:N
        Y(1,I)    =sum(y{1}.*    exp(-a*(X(I)-x{1}).^2));
        Y(2,I)    =sum(y{2}.*    exp(-a*(X(I)-x{2}).^2));
        Y(3,I)    =sum(y{3}.*    exp(-a*(X(I)-x{3}).^2));
        Ynorm(1,I)=sum(ynorm{1}.*exp(-a*(X(I)-x{1}).^2));
        Ynorm(2,I)=sum(ynorm{2}.*exp(-a*(X(I)-x{2}).^2));
        Ynorm(3,I)=sum(ynorm{3}.*exp(-a*(X(I)-x{3}).^2));

      %if I>100
      %y.*exp(-0.001*(X(I)-x).^2)
      %warning('hh')
      %end
    end
    sumY=sumY+Ynorm;


%calculate reper intencities
    modeindind=2;
    for jj=1:numel(reperfreq)
      reperint(jj) = reperint(jj)+ sum(y{modeindind}(find(abs(freqscaled(modeind{modeindind})-reperfreq(jj))<halfwidth)));
      if j~=1
        reperintnorm(jj,j) = reperintnorm(jj,j-1);
      else
        reperintnorm(jj,j) = 0;
      end
      reperintnorm(jj,j) = reperintnorm(jj,j)+ sum(ynorm{modeindind}(find(abs(freqscaled(modeind{modeindind})-reperfreq(jj))<halfwidth)));
    end

    [j reperintnorm(:,j)']
%----------------------------

    Ylim1=-1;
    if moltype==12 %dUrd
      Xlim1=3350;
      Xlim2=3700;
      Ylim2=650; 
    elseif moltype==13 %dThd
      Xlim1=3350;
      Xlim2=3700;
      Ylim2=650; 
    elseif moltype==140 %rAdo
      Xlim1=3250;
      Xlim2=3725;
      Ylim2=0; 
    elseif moltype==540 %rUrd
      Xlim1=3350;
      Xlim2=3700;
      Ylim2=350; 
    else
      Ylim2=0; 
      Xlim1=3350;
      Xlim2=3700;
    end

    if j==1
      Xind=find((X>Xlim1).*(X<Xlim2));
    end

        curplot=0; %current plot number in subplot
    plotsnum=fl_plotcurconf+fl_plotdeconvolsect1+fl_plotdeconvolsect2+1; %number of plots in subplot

%-------------------------------------------1
  %  plot(x,y);
  %  axis([0 1.05*max(x) 0 1.05*max(y)])
    if fl_plotcurconf %plot current conformer spectra
    curplot=curplot+1;

        subplot(plotsnum,1,curplot);
        set(gca,'TickDir','out')
        set(gca,'XMinorTick','on','YMinorTick','on')
        hold on
        if fl_plotgauss
          delete(hgauss1);
          hgauss1 = plot(X,Y(1,:),color(mod(j,numel(color))+1)); %without scaling depending on conformer energy
          set(hgauss1,'LineWidth',2)
        end

        if fl_plotlines
          hline1 = plot([x{1}; x{1}],[zeros(size(ynorm{1})); ynorm{1}],color(mod(j,numel(color))+1));    
        end
        %  axis([0 1.05*max(X_) 0 1.05*max(Y)])
        axis([Xlim1 Xlim2 0 1.05*max(Y(1,Xind))]); %by Y axis unscaled intencity is plotted

  %      xlabel('frequency, cm^-1');
  %      ylabel('IR intensity, KM/Mole');
        title([molname ' ' ms0.desc ' / ' ms0.prop.sdesc ' conformer']);
        box on
    end

%-------------------------------------------2
    if fl_plotdeconvolsect1 %plot first deconvolution section
    curplot=curplot+1;

        subplot(plotsnum,1,curplot);
        set(gca,'TickDir','out')
        set(gca,'XMinorTick','on','YMinorTick','on')
        hold on
        if fl_plotgauss
          delete(hgauss2);
          hgauss2 = plot(X,sumY(2,:),'k'); %without scaling depending on conformer energy
          set(hgauss2,'LineWidth',2)
        end

        if fl_plotlines
          hline2 = plot([x{2}; x{2}],[zeros(size(ynorm{2})); ynorm{2}],color(mod(j,numel(color))+1));    
        end

        if Ylim2
            axis([Xlim1 Xlim2 Ylim1 Ylim2]);
        else
            warning('Cannot find axis scaling data for specified molecule. Defaults used.');
            axis([Xlim1 Xlim2 Ylim1 1.05*max(sumY(1,Xind))]);
        end

    %      xlabel('frequency, cm^-1');
    %      ylabel('IR intensity, KM/Mole');
          title([molname ' summary spectra of ' strcat(labelHval{:}) ' valence mode. ' int2str(j) ' conformer(s) contribution']);
        box on

    end

%-------------------------------------------3
    if fl_plotdeconvolsect2 %plot second deconvolution section
    curplot=curplot+1;
        subplot(plotsnum,1,curplot);
    set(gca,'TickDir','out')
    set(gca,'XMinorTick','on','YMinorTick','on')
        hold on
        if fl_plotgauss
          delete(hgauss3);
          hgauss3 = plot(X,sumY(3,:),'k'); %without scaling depending on conformer energy
          set(hgauss3,'LineWidth',2)
        end

        if fl_plotlines
          hline3 = plot([x{3}; x{3}],[zeros(size(ynorm{3})); ynorm{3}],color(mod(j,numel(color))+1));    
        end

        if Ylim2
            axis([Xlim1 Xlim2 Ylim1 Ylim2]);
        else
            warning('Cannot find axis scaling data for specified molecule. Defaults used.');
            axis([Xlim1 Xlim2 Ylim1 1.05*max(sumY(1,Xind))]);
        end

    %      xlabel('frequency, cm^-1');
    %      ylabel('IR intensity, KM/Mole');
          title([molname ' summary spectra of ' strcat(labelHval2{:}) ' valence mode. ' int2str(j) ' conformer(s) contribution']);
        box on

        %  XTickLabel=get(gca,'XTickLabel');
        %  set(gca,'XTickLabel',XTickLabel(end:-1:1,:));
    end

%-------------------------------------------4 
    curplot=curplot+1;
    subplot(plotsnum,1,curplot);
    set(gca,'TickDir','out')
    set(gca,'XMinorTick','on','YMinorTick','on')

  %  plot(X,sumY,color(mod(j,numel(color))+1));
    hold on
    if fl_plotgauss
      delete(hgauss4);
      hgauss4=plot(X,sumY(1,:),'k');
      set(hgauss4,'LineWidth',2)
    end
    if fl_plotlines
      plot([x{1}; x{1}],[zeros(size(ynorm{1})); ynorm{1}],'k');    
    end



    if Ylim2
        axis([Xlim1 Xlim2 Ylim1 Ylim2]);
    else
        warning('Cannot find axis scaling data for specified molecule. Defaults used.');
        axis([Xlim1 Xlim2 Ylim1 1.05*max(sumY(1,Xind))]);
    end

  %  XTickLabel=get(gca,'XTickLabel');
  %  set(gca,'XTickLabel',XTickLabel(end:-1:1,:));


    xlabel('frequency, cm^-1','FontSize',18);
    ylabel('IR intensity, KM/Mole','FontSize',18);
    box on

  %  confnumber=input('Conformation number?: ');

    title([molname ' IR spectra: ' int2str(j) ' conformer(s) contribution'],'FontSize',18);
  %  pause
    if exist('numconf2plot','var') && any(j==numconf2plot), pause, end;

    j=j+1;

end %cycle over conformers

%print(gcf,'-deps','Filename2.eps')
