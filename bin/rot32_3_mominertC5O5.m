%calculating monent of inertia of H53 in response round C5O5 rotation
%nulibr - frequency of O5H53 libration mode : tbeta
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

moltype=16  %#ok
theory='dftV2'  %#ok
%T=0;
T=298.15;
%----------------------------------

if moltype==9
    nulibr=[271.1304; 203.8675; 272.5326; 254.5251];
    onlyoriginal=0  % process db with only original conformations
    conf2includeind=[99; 4; 103; 2];
    conf2include = [{'EabcA'}; {'AabcA'}; {'EabaA'}; {'AabaA'}];
elseif moltype==12
    nulibr=[279.5929; 210.4989; 278.8513; 254.2816];
    onlyoriginal=0  % process db with only original conformations
    conf2includeind=[54; 7; 52; 4];
    conf2include = [{'EabcA'}; {'AabcA'}; {'EabaA'}; {'AabaA'}]; 
elseif moltype==13
    nulibr=[269,3958; 209,7363; 267,4773; 247,8071; 186,6074; 179,9746 ]; % EabaA-3modes; AabcA-(2-3)modes etc
    onlyoriginal=0  % process db with only original conformations
    conf2includeind=[38; 6; 37; 4; 7; 29];
    conf2include = [{'EabcA'}; {'AabcA'}; {'EabaA'}; {'AabaA'}; {'AabcS'}; {'BabaS'}]; 
elseif moltype==15
    nulibr=[263.7614; 197.4575; 269.9265; 215.8396];
    onlyoriginal=0  % process db with only original conformations
    conf2includeind=[54; 95; 50; 135];
    conf2include = [{'EabcA'}; {'AabcA'}; {'EabaA'}; {'AabaA'}]; 
elseif moltype==16
    nulibr=[259.5264; 186.2395; 262.6992; 181.9675];
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

    iC5  = ms0.ind(find(find(strcmp(pind.labels,'pC5'))==ms0.pind)); 
    iO5 = ms0.ind(find(find(strcmp(pind.labels,'pO5'))==ms0.pind)); 
    iH53 = ms0.ind(find(find(strcmp(pind.labels,'pH53'))==ms0.pind)); 


    j=iH53;
    [C,atomtypenum,IB]=intersect(GLaspec.type,ms0.labels(j));

    m=GLaspec.atommass(atomtypenum)*CC.amu; %[kg]
    IR=m*1E-20*dist2linecart([ms0.x(j) ms0.y(j) ms0.z(j)],[ms0.x(iC5) ms0.y(iC5) ms0.z(iC5)],[ms0.x(iO5) ms0.y(iO5) ms0.z(iO5)])^2; %[kg*m^2]

    ktau = (2*pi*nulibr(indcycle))^2*IR; 


    fi_avg0 = sqrt(CC.k*T/ktau)*360/(2*pi); %[degree]
    warning('off','MATLAB:divideByZero');
    fi_avg = sqrt(CC.h/(8*pi^2*IR*nulibr(indcycle))*coth(CC.h*nulibr(indcycle)/(2*CC.k*T)))*360/(2*pi);
    warning('on','MATLAB:divideByZero');
    

    %                      sm-1                             kcal/rad^2                cm-1
    disp([{ms0.prop.sdesc} nulibr(indcycle)/CC.freqcoef IR ktau/CC.hartree*CC.encoef fi_avg fi_avg0])
    
end

