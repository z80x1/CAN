%
% Version 1.0    
% Last modified  R O Zhurakivsky 2006-03-09
% Created        R O Zhurakivsky 2006-03-09

%EaacS 49
%EaaaS 46
%EabcA 52
%AcbaA 30
%BcbaA 94
%BaacS 39
%BaaaS 37
 
atomsind
T=420

ESP=[
-.83387025610232D+03 
-.83387019646477D+03 
-.83386753566976D+03 
-.83386522241057D+03 
-0.83386482783741D+03
];

ZPE=[
0.223040
0.223043
0.222279 
0.222360 
0.222309
];

dG300=[
.182083	
.182107	
.180311	
.179998	
0.179371
];

dH300=[
.237847	
.237863	
.237595	
.237672	
0.237731
];

dG420=[
0.157170
0.157192
0.154715
0.154161
0.153310
];

fullE=ESP;
dE00all=(fullE-min(fullE))*CC.encoef;

fullE=ESP+ZPE;
dE0all=(fullE-min(fullE))*CC.encoef;

fullE=ESP+dH300;
dH300all=(fullE-min(fullE))*CC.encoef;

fullE=ESP+dG300;
dG300all=(fullE-min(fullE))*CC.encoef;

fullE=ESP+dG420;
dG420all=(fullE-min(fullE))*CC.encoef;

desc=[{'dE00all'} {'dE0all'} {'dH300all'} {'dG300all'} {'dG420all'}]
data=num2cell([dE00all dE0all dH300all dG300all dG420all])


energy=(ESP+dG300-min(ESP+dG300))*CC.hartree;
%energy=(ESP+dG420-min(ESP+dG420))*CC.hartree;
normcoef=1/(sum(exp(-energy/CC.k/T)));
population	=normcoef*exp(-energy/CC.k/T)


Io3h_46=17.96;
Io3h_49=17.84;

energy=(ESP+dG420-min(ESP+dG420))*CC.hartree;
normcoef=1/(sum(exp(-energy(1:2)/CC.k/T)));
pop12=normcoef*exp(-energy(1:2)/CC.k/T)

Io3h_12=Io3h_46*pop12(1) + Io3h_49*pop12(2)  % I=17.90
%--------------------------------

Io5hfree_49=14.3; %87 / EaccS
Io5hfree_46=13.6; %85 / EacaS

energy=(ESP+dG420-min(ESP+dG420))*CC.hartree;
normcoef=1/(sum(exp(-energy(1:2)/CC.k/T)));
pop12=normcoef*exp(-energy(1:2)/CC.k/T)

Io5hfree_12=Io5hfree_49*pop12(1) + Io5hfree_46*pop12(2)  %I=13.96
%--------------------------------

dnu=183  %Iogansen
dH=0.33*sqrt(dnu-40)
dII0=(dH/2.92)^2  %dII0 = 1.8264
%--------------------------------

Io5h_30=40.2
Io5h_94=38.6

energy=(ESP+dG420-min(ESP+dG420))*CC.hartree;
normcoef=1/(sum(exp(-energy(4:5)/CC.k/T)));
pop12=normcoef*exp(-energy(4:5)/CC.k/T)

Io5hfree_12=Io5h_30*pop12(1) + Io5h_94*pop12(2)  %I=


