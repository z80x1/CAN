% take nucleobase are try to form elementary crystal cell
% thymine: a=12.87A, b=6.83A, c=6.70A, beta=105deg, Z=4, P2_1/c
%
% Version 1.0    
% Last modified  R O Zhurakivsky 2007-04-05
% Created        R O Zhurakivsky 2006-?-?

format compact
clear
global pind
global flPlot

atomsind
flPlot=0;

workname='r26001';
odir=[CD.xyzdir filesep workname];
moltype=201;
gtemplname=[workname '_templ.gjf']
fullgtemplname=[CD.templatesdir filesep gtemplname]

if moltype==401
  %Thymine [Oze69] P2_1/c
  cell_a=12.87;
  cell_b=6.83;
  cell_c=6.70;
  cell_beta=105;
  afile = 'r90004_Thymidine.xyz'
elseif moltype==201
  %Cytidine [Bar64] P2_1,2_1,2_1
  cell_a=13.04;
  cell_b=9.49;
  cell_c=3.81;
  afile = 'r90004_Cytidine.xyz'
end

fullafilename=['..\geom\bases' filesep afile]
%aH1_H1 = 10; %indexes of glycosidic bond atoms %H1
%aN1_N1 = 4; %N1

mc0=loadmolxyz(fullafilename);

mc0 = createbondtable(mc0);
mc0 = identmol(mc0,moltype); 
mc0.x = -mc0.x;
mc0.y = -mc0.y;
mc0.z = mc0.z;


cell_a0=0.75;
cell_b0=0.75;
cell_c0=0.75;
cell_da=0.5;
cell_db=0.5;
cell_dc=0.5;

mc1=mc0;
mc1.x = mc1.x  + cell_a0*cell_a;
mc1.y = mc1.y  + cell_b0*cell_b;
mc1.z = mc1.z  + cell_c0*cell_c;
mc2=mc0; %mirror XZ
mc2.x = -mc2.x + (cell_a0-cell_da)*cell_a;
mc2.y = mc2.y  + (cell_b0-cell_db)*cell_b;
mc2.z = -mc2.z + cell_c0*cell_c;
mc3=mc0; %mirror XY
mc3.x = -mc3.x + (cell_a0-cell_da)*cell_a;
mc3.y = -mc3.y + cell_b0*cell_b;
mc3.z = mc3.z  + (cell_c0-cell_dc)*cell_c;
mc4=mc0; %mirror YZ
mc4.x = mc4.x  + cell_a0*cell_a;
mc4.y = -mc4.y + (cell_b0-cell_db)*cell_b;
mc4.z = -mc4.z + (cell_c0-cell_dc)*cell_c;

mcell=mergemol(mc1,mc2);
mcell=mergemol(mcell,mc3);
mcell=mergemol(mcell,mc4);

mcell = createbondtable(mcell);
mcell.desc=[workname '_crystcell'];

if 0
    mcell2 = mcell;
    mcell2.x = mcell2.x +cell_a;
    mcell3 = mcell;
    mcell3.y = mcell3.y +cell_b;
    mcell4 = mcell;
    mcell4.z = mcell4.z +cell_c;
    mcell5 = mcell;
    mcell5.x = mcell5.x +cell_a;
    mcell5.y = mcell5.y +cell_b;
    mcell6 = mcell;
    mcell6.y = mcell6.y +cell_b;
    mcell6.z = mcell6.z +cell_c;
    mcell7 = mcell;
    mcell7.x = mcell7.x +cell_a;
    mcell7.z = mcell7.z +cell_c;
    mcell8 = mcell;
    mcell8.x = mcell8.x +cell_a;
    mcell8.y = mcell8.y +cell_b;
    mcell8.z = mcell8.z +cell_c;

    mcellf=mergemol(mcell,mcell2);
    mcellf=mergemol(mcellf,mcell3);
    mcellf=mergemol(mcellf,mcell5);
    %mcellf=mergemol(mcellf,mcell4);
    %mcellf=mergemol(mcellf,mcell6);
    %mcellf=mergemol(mcellf,mcell7);
    %mcellf=mergemol(mcellf,mcell8);
    mcellf = createbondtable(mcellf);
end

if exist(odir)~=7
    mkdir(odir);
end
savemol(odir,mcell,0);

order=[];
savemolgs(odir,mcell,3,order,fullgtemplname); %Gaussian with XYZ

plotmol(mcell)
%plotmol(mcellf)


