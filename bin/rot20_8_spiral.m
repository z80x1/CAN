%make spiral from one sugar residue

%
% Version 1.0    
% Last modified  R O Zhurakivsky 2009-02-09
% Created        R O Zhurakivsky 2009-02-09

%molA - first molecule
%molZ - resulting molecule

%handle=figure('Color','White');
%hold on


format compact
clear
global pind
global flPlot

flPlot=1;

odir=['../xyz' filesep 'r90201']; %output directory
molAfname = 'r902_C5PH9O6.xyz' %mol_B
fullmolAfname=['../templ' filesep molAfname]

iA_begin=uint16([18 14 9]); %
iA_end=uint16([5 1 19]); %

junk_val1 = 119.5/180*pi; % iA_end(2) - juncion - iA_begin(2)
junk_tor1 = 149.769/180*pi;  % iA_end(1) - iA_end(2) - juncion - iA_begin(2)

junk_tor2 = -87.031/180*pi;  % iA_end(2) - juncion - iA_begin(2) - iA_begin(3)



    molA0=loadmolxyz(fullmolAfname);
    indmove=0;

    iA_begin2Bis= iA_begin;
    iA_endBis = iA_end;
    molZ = molA0;

for i=1:10


plotmol(molA0,'r',0);

    iA_beginBis=iA_begin2Bis;
    iA_begin2Bis=iA_beginBis+molA0.atomnum;

    if i~=1
        iA_endBis=iA_endBis+molA0.atomnum;
    end

    %movement vector from begin to end atoms
    molB0=molA0;
    vect_BE = createvect(molB0, iA_begin(1), molZ, iA_endBis(3))
    molB0.x = molB0.x + vect_BE(1);
    molB0.y = molB0.y + vect_BE(2);
    molB0.z = molB0.z + vect_BE(3);

    molZ0=mergemol(molZ,molB0);

    fi = p2pangle(createplane(molZ0,iA_endBis(1),molZ0,iA_endBis(2),molZ0,iA_endBis(3)),...
                  createplane(molZ0,iA_endBis(2),molZ0,iA_begin2Bis(1),molZ0,iA_begin2Bis(2)));
    fi2 = p2pangle(createplane(molZ0,iA_endBis(2),molZ0,iA_begin2Bis(1),molZ0,iA_begin2Bis(2)),...
                  createplane(molZ0,iA_begin2Bis(1),molZ0,iA_begin2Bis(2),molZ0,iA_begin2Bis(3)));
%fi*180/pi
%fi2*180/pi


    vect_B = createvect(molB0, iA_end(2), molB0, iA_end(3));
    iB0_center = iA_begin(1);
    molB1 = rotline3(vect_B,molB0,junk_tor1-fi,iB0_center);
%plotmol(molB0,'b',0);
plotmol(molB1,'k',0);

    vect_B = createvect(molB1, iA_begin(1), molB1, iA_begin(2));
    iB0_center = iA_begin(1);
    molB2 = rotline3(vect_B,molB1,junk_tor2-fi2,iB0_center);
plotmol(molB2,'r',0);

%    molB3 = delatom( molB2, iA_begin(1) );
    molB3=molB2;

    molZ=mergemol(molZ,molB3);

end    

    molZ.desc=['spiral'];
    molZ = createbondtable(molZ);

%plotmol(molZ);

    if exist(odir)~=7
        mkdir(odir);
    end
    savemol(odir,molZ,0);



