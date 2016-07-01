%calculating monent of inertia of sugar and base in response to glicosidic bond 
%nulibr - frequency of C5HHO5H53 and other part of nucleoside round C4C5 libration mode : tgamma
%
% Version 1.0    
% Last modified  R O Zhurakivsky 2008-11-24
% Created        R O Zhurakivsky 2006-?-?

format compact
clear
global pind
global flPlot

pindsdef
atomsind

%----------------------------------
flPlot=0;

moltype=13  %#ok
theory='dftV2'  %#ok
%T=0;
T=298.15;
%----------------------------------

if moltype==9
    nulibr=[123.9; 123.8; 125.8; 120.1];
    onlyoriginal=0  % process db with only original conformations
    conf2includeind=[99; 4; 103; 2];
    conf2include = [{'EabcA'}; {'AabcA'}; {'EabaA'}; {'AabaA'}];
elseif moltype==12
    nulibr=[136.6; 125.3; 136.8; 122.9];
    onlyoriginal=0  % process db with only original conformations
    conf2includeind=[54; 7; 52; 4];
    conf2include = [{'EabcA'}; {'AabcA'}; {'EabaA'}; {'AabaA'}]; 
elseif moltype==13
    nulibr=[142.36; 118.36; 136.96; 111.94]; %no modes!! (it's difficult to distiguish librational modes around C4C5 bond)
    onlyoriginal=0  % process db with only original conformations
    conf2includeind=[38; 6; 37; 4];
    conf2include = [{'EabcA'}; {'AabcA'}; {'EabaA'}; {'AabaA'}]; 
elseif moltype==15
    nulibr=[135.7; 112.8; 134.8; 112.6];
    onlyoriginal=0  % process db with only original conformations
    conf2includeind=[54; 95; 50; 135];
    conf2include = [{'EabcA'}; {'AabcA'}; {'EabaA'}; {'AabaA'}]; 
elseif moltype==16
    nulibr=[134.3; 137.2; 132.8; 134.2];
    onlyoriginal=0  % process db with only original conformations % EabaA-2modes AabaA-2modes
    conf2includeind=[55; 102; 51; 100];
    conf2include = [{'EabcA'}; {'AabcA'}; {'EabaA'}; {'AabaA'}]; 
else
    error('Dont know anything about conf2include & nulibr');
end
nulibr=nulibr*CC.freqcoefB3LYP631Gdp*CC.freqcoef; %[Hz=1/c]


if ~strcmp(theory,'dft')
  theorystr = ['_' theory];
else                                                                                    
  theorystr = '';
end
workdbname=[CD.dbdir filesep 'r' int2str(moltype) '_g' theorystr];
if onlyoriginal
    templ='_or';
    workdbname = [workdbname templ];
end
workdbname=[workdbname '.mat']  %#ok

load(workdbname,'workdb')

recnum=numel(workdb);
for i=1:recnum
    sdesc(i) = {workdb(i).prop.sdesc};
end

if isempty(conf2includeind)
    ind=[];
    for i=1:size(conf2include,1)
        ind(end+1)=strcmpcellar(sdesc,conf2include(i,:));
    end
else
    ind=conf2includeind;
end

disp(['ms0.prop.sdesc  ' 'nulibr  ' 'IR  ' 'ktau  ' 'fi_avg  ' 'fi_avg0  ' 'fi_avg00  ' 'nu_Isugar_inf  '])
indcycle=0;
ind=reshape(ind,1,numel(ind));
for i=ind %sorting conformers
    indcycle=indcycle+1;
    ms0=workdb(i);

    iC4  = ms0.ind(find(find(strcmp(pind.labels,'pC4'))==ms0.pind)); 
    iC5  = ms0.ind(find(find(strcmp(pind.labels,'pC5'))==ms0.pind)); 

    Ires1=0;
    Ires2=0;

    indres1=[];
    indres2=[];
    for j=1:ms0.atomnum

        if strcmp(pind.labels{ms0.pind(j)},'pH51')
            indres1(1)=j;
        elseif strcmp(pind.labels{ms0.pind(j)},'pH52')
            indres1(2)=j;
        elseif strcmp(pind.labels{ms0.pind(j)},'pH53')
            indres1(3)=j;
        elseif strcmp(pind.labels{ms0.pind(j)},'pO5')
            indres1(4)=j;
        else
            indres2(end+1)=j;
        end
    end
    if sum(indres1==0)
        error('not all of indexes found!');
    end
    
    for j=indres1
        [C,atomtypenum,IB]=intersect(GLaspec.type,ms0.labels(j));

        m=GLaspec.atommass(atomtypenum)*CC.amu; %[kg]
        Ires1=Ires1+m*1E-20*dist2linecart([ms0.x(j) ms0.y(j) ms0.z(j)],[ms0.x(iC4) ms0.y(iC4) ms0.z(iC4)],[ms0.x(iC5) ms0.y(iC5) ms0.z(iC5)])^2; %[kg*m^2]
    end
    for j=indres2
        [C,atomtypenum,IB]=intersect(GLaspec.type,ms0.labels(j));

        m=GLaspec.atommass(atomtypenum)*CC.amu;
        Ires2=Ires2+m*1E-20*dist2linecart([ms0.x(j) ms0.y(j) ms0.z(j)],[ms0.x(iC4) ms0.y(iC4) ms0.z(iC4)],[ms0.x(iC5) ms0.y(iC5) ms0.z(iC5)])^2;
    end

    IR=1/(1/Ires1+1/Ires2); %[kg*m^2] %total moment of inertia
    ktau = (2*pi*nulibr(indcycle))^2*IR; %[J/rad^2] rotary spring force

    ktau0 = (2*pi*nulibr(indcycle))^2*Ires2;  %[J/rad^2] rotary spring force when Ires1=infinity
    fi_avg00 = sqrt(CC.k*T/ktau0)*360/(2*pi); %[degree] fi amplitude in k_tau*fi^2/2=kT/2 assumption


    fi_avg0 = sqrt(CC.k*T/ktau)*360/(2*pi); %[degree] fi amplitude in k_tau*fi^2/2=kT/2 assumption when Ires1=infinity
    % <fi^2>=h/(2*pi)/(2*IR*omega)*cth(h*nu/2kT)
    warning('off','MATLAB:divideByZero');
    fi_avg = sqrt(CC.h/(8*pi^2*IR*nulibr(indcycle))*coth(CC.h*nulibr(indcycle)/(2*CC.k*T)))*360/(2*pi);
    warning('on','MATLAB:divideByZero');
    
    nu_Ires1_inf=1/2/pi*sqrt(ktau/Ires2); %[1/s] mode frequency in case of  Ires1=infinity

    %                       sm-1                                           kcal/rad^2                                        cm-1
    disp([{ms0.prop.sdesc} nulibr(indcycle)/CC.freqcoef IR ktau/CC.hartree*CC.encoef fi_avg fi_avg0 fi_avg00 nu_Ires1_inf/CC.freqcoef])
    
end

                                                                                          
