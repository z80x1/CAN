%calculating monent of inertia of sugar and base in response to glicosidic bond 
%nulibr - frequency of sugar and base round C1N1 libration mode
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

moltype=15  %#ok
theory='dftV2'  %#ok
T=0;
%T=298.15;
%----------------------------------

if moltype==9
    nulibr=[49,46; 22,03; 52,36; 33,58; 44,52; 44,52; 30,29; ];
%    nulibr=[3237; 3235; 3236; 3228; 3233; 3233];
    onlyoriginal=0  % process db with only original conformations
    conf2includeind=[4; 99; 2; 103; 66; 47; 139];
    conf2include = [{'AabcA'}; {'EabcA'}; {'AabaA'}; {'EabaA'}; {'AabcS'}; {'BabaS'}; {'CabcA'}];
elseif moltype==12
    nulibr=[41.15; 45.3; 44.81; 44.72];
    onlyoriginal=0  % process db with only original conformations
    conf2includeind=[54; 7; 52; 4];
    conf2include = [{'EabcA'}; {'AabcA'}; {'EabaA'}; {'AabaA'}]; 
elseif moltype==13
    nulibr=[34,68; 37,39; 32,59; 38,40; 42,03; 39,74; 31,63; ];
%    nulibr=[3230; 3214; 3229; 3212; 3209; 3210];
    onlyoriginal=0  % process db with only original conformations
    conf2includeind=[6; 38; 4; 37; 7; 29; 125];
    conf2include = [{'AabcA'}; {'EabcA'}; {'AabaA'}; {'EabaA'}; {'AabcS'}; {'BabaS'}; {'CabcA'}]; 
elseif moltype==15
    nulibr=[31,59; 32,39; 15,56; 34,14; 25,86];
%    nulibr=[3266; 3260; 3265; 3267; 3264];
    onlyoriginal=0  % process db with only original conformations
    conf2includeind=[54; 142; 95; 50; 135];
    conf2include = [{'EabcA'}; {'CabcA'}; {'AabcA'}; {'EabaA'}; {'AabaA'}]; 
elseif moltype==16
    nulibr=[30,70; 26,43; 19,58; 30,85; 21,80];
%    nulibr=[3271; 3265; 3270; 3272; 3271];
    onlyoriginal=0  % process db with only original conformations
    conf2includeind=[55; 141; 102; 51; 100];
    conf2include = [{'EabcA'}; {'CabcA'}; {'AabcA'}; {'EabaA'}; {'AabaA'}]; 
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

    iC1  = ms0.ind(find(find(strcmp(pind.labels,'pC1'))==ms0.pind)); 
    if any(moltype==[11,9,12,13,14,17,18,19,21]) ||...
        (moltype>=200 && moltype<=299) || (moltype>=400 && moltype<=499) || (moltype>=500 && moltype<=599) %pirimidines
      i2 = ms0.ind(find(find(strcmp(pind.labels,'bN1'))==ms0.pind)); 
    elseif any(moltype==[15,16]) || (moltype>=100 && moltype<=199) || (moltype>=300 && moltype<=399) %purines
      i2 = ms0.ind(find(find(strcmp(pind.labels,'bN9'))==ms0.pind)); 
    else
      error('Cannot locate atoms of glycosidic bond');
    end

    Isugar=0;
    Ibase=0;

    indsugar=[];
    indbase=[];
    for j=1:ms0.atomnum
        if pind.labels{ms0.pind(j)}(1)=='p'
            indsugar(end+1)=j;
        end
        if pind.labels{ms0.pind(j)}(1)=='b'
            indbase(end+1)=j;
        end
    end
       
    for j=indsugar
        [C,atomtypenum,IB]=intersect(GLaspec.type,ms0.labels(j));

        m=GLaspec.atommass(atomtypenum)*CC.amu; %[kg]
        Isugar=Isugar+m*1E-20*dist2linecart([ms0.x(j) ms0.y(j) ms0.z(j)],[ms0.x(iC1) ms0.y(iC1) ms0.z(iC1)],[ms0.x(i2) ms0.y(i2) ms0.z(i2)])^2; %[kg*m^2]
%[j pind.labels(ms0.pind(j)) ms0.x(j) ms0.y(j) ms0.z(j) dist2linecart([ms0.x(j) ms0.y(j) ms0.z(j)],[ms0.x(iC1) ms0.y(iC1) ms0.z(iC1)],[ms0.x(i2) ms0.y(i2) ms0.z(i2)])^2 ]
%disp(sprintf('%s %0.5g %0.5g %0.5g',ms0.labels{j}, ms0.x(j), ms0.y(j), ms0.z(j)))
    end
    for j=indbase
        [C,atomtypenum,IB]=intersect(GLaspec.type,ms0.labels(j));

        m=GLaspec.atommass(atomtypenum)*CC.amu;
        Ibase=Ibase+m*1E-20*dist2linecart([ms0.x(j) ms0.y(j) ms0.z(j)],[ms0.x(iC1) ms0.y(iC1) ms0.z(iC1)],[ms0.x(i2) ms0.y(i2) ms0.z(i2)])^2;
%[j pind.labels(ms0.pind(j)) ms0.x(j) ms0.y(j) ms0.z(j) dist2linecart([ms0.x(j) ms0.y(j) ms0.z(j)],[ms0.x(iC1) ms0.y(iC1) ms0.z(iC1)],[ms0.x(i2) ms0.y(i2) ms0.z(i2)])^2 ]
%disp(sprintf('%s %0.5g %0.5g %0.5g',ms0.labels{j}, ms0.x(j), ms0.y(j), ms0.z(j)))
    end


    IR=1/(1/Isugar+1/Ibase); %[kg*m^2] %total moment of inertia
    ktau = (2*pi*nulibr(indcycle))^2*IR; %[N*m/rad^2] rotary spring force

    ktau0 = (2*pi*nulibr(indcycle))^2*Ibase;  %[N*m/rad^2] rotary spring force when Isugar=infinity
    fi_avg00 = sqrt(CC.k*T/ktau0)*360/(2*pi); %[degree] fi amplitude in k_tau*fi^2/2=kT/2 assumption when Isugar=infinity


    fi_avg0 = sqrt(CC.k*T/ktau)*360/(2*pi); %[degree] fi amplitude in k_tau*fi^2/2=kT/2 assumption 
    % <fi^2>=h/(2*pi)/(2*IR*omega)*cth(h*nu/2kT)
    warning('off','MATLAB:divideByZero');
    fi_avg = sqrt(CC.h/(8*pi^2*IR*nulibr(indcycle))*coth(CC.h*nulibr(indcycle)/(2*CC.k*T)))*360/(2*pi);
    warning('on','MATLAB:divideByZero');

    nu_Isugar_inf=1/2/pi*sqrt(ktau/Ibase); %[1/s] mode frequency in case of  Isugar=infinity

    %                       sm-1                                           kcal/rad^2                                        cm-1
    disp([{ms0.prop.sdesc} nulibr(indcycle)/CC.freqcoef IR ktau/CC.hartree*CC.encoef fi_avg fi_avg0 fi_avg00 nu_Isugar_inf/CC.freqcoef])

end

                                                                                          