function d = bondlen(a1,a2)
%returns typical bond length between atoms of types a1 & a2
%
% Version 1.0    
% Last modified  R O Zhurakivsky 2008-04-22 
% Created        R O Zhurakivsky 2005-09-?

%[a1,a2]=sort([a1,a2]);

%values from Seanger
%CC 1.537
%CN 1.475
%CO 1.426
%CF 1.379
%CP 1.841
%PO 1.71

%values of covalent radia from "Таблица Менделеева 2.0" [http://www.kornsoft.narod.ru/progi/mendel.htm]
%value for Na is reduced from 1.54 to 1 for preventing bonds creation
corad.labels=[{'C'},{'N'},{'O'},{'P'},{'F'},{'H'},{'Na'},{'Br'},{'S'}];
%corad.value =[0.77  0.75  0.73  1.06  0.72  0.32  1.00];
corad.value =[0.77  0.75  0.73  1.06  0.72  0.40  1.00    1.14  1.02]; %071018 value for H changed from 0.32 to 0.40



d = 1.1*(corad.value(strcmpcellar(corad.labels,a1))+corad.value(strcmpcellar(corad.labels,a2)));
