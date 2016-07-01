%rot27: calculate vibrational, rotational energies and entropy
%not_completed!!
%
% Version 1.0    
% Last modified  R O Zhurakivsky 2008-05-19
% Created        R O Zhurakivsky 2006-?-?

clear
format compact

atomsind
%T=298.15
T=1E-6
p=101325; %presure in Pascals
n=1; %number of mols

%IA=1532.85391*CC.amu*CC.l^2;
%IB=5752.71767*CC.amu*CC.l^2;
%IC=6304.30599*CC.amu*CC.l^2;



    We=1; %electron state degeneracy
    Se = n*CC.R*log(We);
    Htrans = 3/2*CC.R*T;
    Hrot = 3/2*CC.R*T;


workdbname=[CD.dbdir filesep 'r12_g_test.mat']
load(workdbname,'workdb')

    M = 0; % molecule mass (kg)
    for jj=1:numel(GLaspec.type)
	M = M+GLaspec.atommass(jj)*sum(workdb(1).labels==GLaspec.type(jj));
    end	
    M = M*CC.amu;

    Strans = n*CC.R*(1.5+log((pi*M*CC.k*T)^1.5*(n*CC.R*T/p)));

recnum=numel(workdb);
for i=1:recnum
  if workdb(i).new=='Y'

    ms0=workdb(i);


    sdesc(i,:)=ms0.prop.sdesc;

    ZPEg(i)=ms0.gaussian.ZPE;
    TECg(i)=ms0.gaussian.TEC;
    EECg(i)=ms0.gaussian.EEC;
    GECg(i)=ms0.gaussian.GEC;
    freqs=ms0.freq.freq(find(ms0.freq.freq>0))*CC.freqcoef;
%    E0 = ms0.gaussian.MP2_6311__Gdp*CC.hartree;
    E0=0;

    imoments = ms0.gaussian.imoments;
    VA=CC.h^2/(8*pi*imoments(1)*CC.k*T);
    VB=CC.h^2/(8*pi*imoments(2)*CC.k*T);
    VC=CC.h^2/(8*pi*imoments(3)*CC.k*T);


    ZPE(i)=CC.h/2*sum(freqs)/CC.hartree*CC.encoef; %kcal/mol


%    EvibT(i) = CC.k*T*log(prod(1-exp(-CC.h*freqs/(CC.k*T)))) /CC.hartree*CC.encoef; %kcal/mol

%    Erot(i)=-CC.k*T*log(pi^(1/3)*((8*pi^2*CC.k*T)/CC.h^2)^1.5*(IA*IB*IC)^1.5)/CC.hartree*CC.encoef;

 
    Ua=CC.h*freqs/(CC.k*T);

    Srot = n*CC.R*(1.5+log(sqrt(pi*VA*VB*VC)/2));
    Svib = n*CC.R*sum((1./(Ua.*exp(Ua)-1))-log(1-exp(-Ua)));
    Sfull(i) = Strans + Srot + Svib + Se - n*CC.R*(log(n*CC.NA)-1);


    Hvib = CC.NA*CC.h*sum(freqs./(exp(Ua)-1));
    Hfull(i) = Htrans + Hrot + Hvib + E0 + CC.R*T;

    G(i)=Hfull(i)-T*Sfull(i);
                     
  end
end

sdesc
ZPEg
TECg
EECg
GECg
ZPE
%EvibT
%Erot
