%calculating monent of inertia of H32 in response round C3O3 rotation
%nulibr - frequency of O3H32 libration mode : tepsilon
%
% Version 1.0    
% Last modified  R O Zhurakivsky 2008-05-19
% Created        R O Zhurakivsky 2006-?-?

format compact
clear
global pind
global flPlot

pindsdef
atomsind

%----------------------------------
flPlot=0;

moltype=16 %#ok
theory='dftV2'  %#ok
%T=0;
T=298.15;
%----------------------------------

if moltype==9
    nulibr=[312.8153; 284.4239; 301.9338; 234.1317];
    onlyoriginal=0  % process db with only original conformations
    conf2includeind=[99; 4; 103; 2];
    conf2include = [{'EabcA'}; {'AabcA'}; {'EabaA'}; {'AabaA'}];
elseif moltype==12
    nulibr=[317.477; 269.476; 305.6015; 235.8666];
    onlyoriginal=0  % process db with only original conformations
    conf2includeind=[54; 7; 52; 4];
    conf2include = [{'EabcA'}; {'AabcA'}; {'EabaA'}; {'AabaA'}]; 
elseif moltype==13
    nulibr=[318,4172; 267,6585; 302,6599; 198,3927; 253,1724; 124,4813];
    onlyoriginal=0  % process db with only original conformations
    conf2includeind=[38; 6; 37; 4; 7; 29];
    conf2include = [{'EabcA'}; {'AabcA'}; {'EabaA'}; {'AabaA'}; {'AabcS'}; {'BabaS'}]; 
elseif moltype==15
    nulibr=[312.1877; 273.3995; 301.3079; 240.7294];
    onlyoriginal=0  % process db with only original conformations
    conf2includeind=[54; 95; 50; 135];
    conf2include = [{'EabcA'}; {'AabcA'}; {'EabaA'}; {'AabaA'}]; 
elseif moltype==16
    nulibr=[311.7686; 255.9778; 306.0207; 242.4575];
    onlyoriginal=0  % process db with only original conformations
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

disp(['ms0.prop.sdesc  ' 'nulibr  ' 'IR  ' 'ktau  ' 'fi_avg  ' 'fi_avg0  ' ])
indcycle=0;
ind=reshape(ind,1,numel(ind));
for i=ind %sorting conformers
    indcycle=indcycle+1;
    ms0=workdb(i);

    iC3  = ms0.ind(find(find(strcmp(pind.labels,'pC3'))==ms0.pind)); 
    iO3 = ms0.ind(find(find(strcmp(pind.labels,'pO3'))==ms0.pind)); 
    iH32 = ms0.ind(find(find(strcmp(pind.labels,'pH32'))==ms0.pind)); 


    j=iH32;
    [C,atomtypenum,IB]=intersect(GLaspec.type,ms0.labels(j));

    m=GLaspec.atommass(atomtypenum)*CC.amu; %[kg]
    IR=m*1E-20*dist2linecart([ms0.x(j) ms0.y(j) ms0.z(j)],[ms0.x(iC3) ms0.y(iC3) ms0.z(iC3)],[ms0.x(iO3) ms0.y(iO3) ms0.z(iO3)])^2; %[kg*m^2]

    ktau = (2*pi*nulibr(indcycle))^2*IR; 


    fi_avg0 = sqrt(CC.k*T/ktau)*360/(2*pi); %[degree]
    warning('off','MATLAB:divideByZero');
    fi_avg = sqrt(CC.h/(8*pi^2*IR*nulibr(indcycle))*coth(CC.h*nulibr(indcycle)/(2*CC.k*T)))*360/(2*pi);
    warning('on','MATLAB:divideByZero');
    

    %                      sm-1                             kcal/rad^2                cm-1
    disp([{ms0.prop.sdesc} nulibr(indcycle)/CC.freqcoef IR ktau/CC.hartree*CC.encoef fi_avg fi_avg0])
    
end

