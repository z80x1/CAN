%
% Version 1.0    
% Last modified  R O Zhurakivsky 2009-01-09 
% Created        R O Zhurakivsky 2005-09-?

format compact
global GLsdesc GLsdescr GLaspec CG CC CD CH VDV magictorsion

GLsdesc=['a','b','c','d','e','f','g','h','i'];
GLsdescr=['A','B','C','D','E','F','G','H','I','J'];

GLaspec.type=    [{'H'};{'C'};{'N'};{'O'};{'F'};{'P'};{'Na'};{'Br'};{'S'};{'Cl'}];
GLaspec.weight=  [1;6;7;8;9;15;11;35;16;17];
GLaspec.atommass=[1.00783; 12; 14.00307; 15.99491; 18.99840; 30.97376; 22.98977; 79.904; 32.066; 35.4527]; %amu (from Gaussian)
%Br atom mass 79.904 need to be checked!!!

%091119 values from Ковалентный радиус: http://ru.wikipedia.org/wiki/Категория:Химические_элементы
GLaspec.corad =  [0.40; 0.77; 0.75; 0.73; 0.72; 1.06; 1.54;  1.14;  1.02; 0.99; ]; 
%values of covalent radia from "Таблица Менделеева 2.0" [http://www.kornsoft.narod.ru/progi/mendel.htm]
%value for Na is reduced from 1.54 to 1 for preventing bonds creation
%071018 value for H changed from 0.32 to 0.40

%bond lengths values from Seanger
%CC 1.537
%CN 1.475
%CO 1.426
%CF 1.379
%CP 1.841
%PO 1.71



%critical values for AH...B bonds
%CG.lcritABdist=3.8;
CG.lcritABdist=3.96; %090604 вибірка по 2583 конформерам 26 молекул
%CG.lcritHBdist=2.72; %060608 Van_der_Vaals radii (Bondi,1964 - Seanger)
CG.lcritHBdist=3.18; %090604 вибірка по 2583 конформерам 26 молекул
%CG.lcritAHBang=90;
CG.lcritAHBang=98.5; %090604 вибірка по 2583 конформерам 26 молекул

CC.encoef = 627.5095; %a.e. to kcal/mol energy conversion coefficient (from Gaussian)
CC.ZPEcoefHF631Gd=0.8929;  %coefficient to convert ZPE energy 

CC.ZPEcoefB3LYP631Gd=0.9804;  %coefficient to convert ZPE energy [Foresman Exploring chemistry with electronic struture methods]
CC.freqcoefB3LYP631Gd=0.9613;  %coefficient to convert mode's frequinces

CC.ZPEcoefB3LYP631Gdp=0.9804;  %coefficient to convert ZPE energy 
CC.freqcoefB3LYP631Gdp=0.9608;  %coefficient to convert mode's frequinces [http://cccbdb.nist.gov/vsf.asp 2010-0407]

CC.freqcoefMP2FC631Gdp=0.9365;

CD.templatesdir=['..' filesep 'templ'];
CD.bindir='.';
CD.dbdir=['..' filesep 'db'];
CD.xlsdir=['..' filesep 'xls'];
CD.datadir=['..' filesep '..' filesep 'CAN2' filesep 'data'];
CD.xyzdir=['..' filesep 'xyz'];
CD.geomdir=['..' filesep 'geom'];
CD.AIM=['..' filesep 'AIMres'];

CC.h=6.6254E-34; %Plank constant [J*s]
CC.k=1.380658E-23; %Boltzman contsant [J/K]
CC.e=1.602E-19; %charge of electron [C]
CC.amu=1.66035E-27; %amu [kg]
CC.me=5.4876E-4*CC.amu; %mass of electron [kg]
CC.R=8.31436; %Gas contant [J/K/mol]
CC.NA=6.0247E23; %Avogadro constant [1/mol]
CC.c=2.997928E8; %speed of light [m/sec]

CC.l=5.29173E-11; %radius of the first Bohr's orbit [m]
CC.rydberg = 3.288E15; %rydberg constant [sec^-1]  [hartree=2*CC.rydrerg*CC.h]


CC.hartree=4.3597482E-18; %[Joules]
CC.freqcoef=CC.c*100; %cm/s - speed of light   
%(kT = 208.53cm-1 = 0.0259eV = 0.5962kcal/mol = 9.5004E-4 Hartree , T=300K)

%chemical constants
CH.typOHNdist=1.7; %typical OH...N distance between NH3 molecule and nucleoside


%Wan-der-vaals radii (13.	ћ.‚. ‡ҐдЁа®ў, Џ.Њ. ‡®аЄЁ© ‚ ­-¤Ґа-‚  «мб®ўл а ¤Ёгбл Ё Ёе ЇаЁ¬Ґ­Ґ­ЁҐ ў еЁ¬ЁЁ // “бЇҐеЁ еЁ¬ЁЁ. - 1989. - LVIII, ўлЇ.5. - ‘.713-746.):
%C = 1.71A
%H = 1.16A
%O = 1.29A
%N = 1.50A
%Cl = 1.90A
%S = 1.84A

%Wan-der-vaals radii
%A. Bondi, J. Phys. Chem., 1964, 68, 441. ??
%http://www.webelements.com/webelements/scholar/elements/
%C = 1.70A
%H = 1.20A
%O = 1.52A
%N = 1.55A
%Cl = 1.75A
%S = 1.80A

%A. Bondi, J. Phys. Chem., 1964, 68, 441. ??
%HO=2.72; HN=2.75; HC=2.90; HH=2.40; 

VDV.Bondi.C=1.70;
VDV.Bondi.H=1.20;
VDV.Bondi.O=1.52;
VDV.Bondi.N=1.55;
VDV.Bondi.Cl=1.75;
VDV.Bondi.S=1.80;

VDV.Zorkiy.C=1.71;
VDV.Zorkiy.H=1.16;
VDV.Zorkiy.O=1.29;
VDV.Zorkiy.N=1.50;
VDV.Zorkiy.Cl=1.90;
VDV.Zorkiy.S=1.84;


%'new' attribute in molecule data structure legend:
%Y	new conformation
%N	similar to one exists
%B	bad: beta-ribose conformation or something else
%P	similar to one in the previous conformation family
%X	similar to one in the next conformation family
%T	transition state
%F	with constrained values

magictorsion(1,:)=[ 6  1 20 21]; %tchi sorted by pind for purines
magictorsion(2,:)=[ 6  1 23 41]; %tchi sorted by pind for pirimidined
magictorsion(3,:)=[19  9  5  4]; %tbeta sorted by pind
magictorsion(4,:)=[ 5  4  3  8]; %tgamma sorted by pind
magictorsion(5,:)=[ 3  4  8 15]; %tepsilon sorted by pind

%Conformation codes
%Y	new conformation
%N	similar to one exists
%B	bad: beta-ribose conformation or something else
%P	similar to one in the previous conformation family
%X	similar to one in the next conformation family
%T	transition state
%D	constrained



