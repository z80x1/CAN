%
% Version 1.0    
% Last modified  R O Zhurakivsky 2007-10-22
% Created        R O Zhurakivsky 2006-?-?

format compact
clear
global pind
global flPlot

atomsind

%----------------------------------
flPlot=0;

%workdir='D:\_diplom\work\Nata\bases\';
%fname='r46002_Thymine_gerdil.xyz';  moltype=460
%fname='r26001_Cytosine_Barker.xyz'; moltype=260
%fname='r26002_Cytosine_Jeffrey.xyz'; moltype=261

%dUrd Rahman 1972
%workdir='D:\_diplom\work\0611-0706.dUrd_eng\070803_crystal\';
%fname='dUrd_crdata_rahman.CCgeom'; moltype=12
%workdir='D:\_diplom\CAN\xyz\rc34\';
%fname='crdata_rahman_withH.xyz'

%dAdo Klooster 1991
%workdir='D:\_diplom\docs\articles\_crystalstructure';
%fname='rAdo_crdata_klooster.CCgeom'; moltype=140

%dAdo Klooster 1991
%workdir='D:\_diplom\docs\articles\_crystalstructure';
%fname='dCyd_crdata_subramanian.CCgeom'; moltype=9

%dGuo Haschmeyer 1965
workdir='D:\_diplom\docs\articles\_crystalstructure\dGuo';
%fname='dGuo_crdata_Haschmeyer.CCgeom'; moltype=16
fname='Br5dCyd-dGuo_crdata_Haschmeyer.CCgeom'; moltype=16

fl_repcells = 0; %build 8 cells in 3D
fl_writefile=1
%----------------------------------

dlm=strfind(fname,'.');
fnameshort = fname(1:(dlm-1));


mc0=loadmolxyz([workdir filesep fname],fnameshort); %M2


if moltype==460 %r46002_Thymine_gerdil.xyz
    cell_a=6.077;
    cell_b=27.862;
    cell_c=3.816;
    cell_beta=94+19/60;

    CS1 = [1/2 0 1/2]; %center of symmetry
    DPy= 1/4; %glide plane  y=1/4

    DAx=1/2; %screw diad axis [1/2 y 3/4]
    DAz=3/4;


    mc1 = mc0; %M3
    mc1.x = -(mc1.x-CS1(1)) + CS1(1);
    mc1.y = -(mc1.y-CS1(2)) + CS1(2);
    mc1.z = -(mc1.z-CS1(3)) + CS1(3);

    mc2 = mc0; %M5
    mc2.x = mc2.x;
    mc2.y = -(mc2.y-DPy) + DPy ; %
    mc2.z = mc2.z + 1/2;

    mc3 = mc2; %
    mc3.x = -(mc3.x-CS1(1)) + CS1(1);
    mc3.y = -(mc3.y-CS1(2)) + CS1(2);
    mc3.z = -(mc3.z-CS1(3)) + CS1(3);

    %mc4 = mc1; %M5
    %mc4.x = -(mc4.x-DAx)+DAx;
    %mc4.y = mc4.y + 1/2;
    %mc4.z = -(mc4.z-DAz)+DAz;


    mcell=mergemol(mc0,mc1);
    mcell=mergemol(mcell,mc2);
    mcell=mergemol(mcell,mc3);
    %mcell=mergemol(mcell,mc4);


    mcell.x=mcell.x*cell_a;  %denormalizing coordinates
    mcell.y=mcell.y*cell_b;
    mcell.z=mcell.z*cell_c;

    mcell.x = mcell.x + mcell.z.*cosd(cell_beta);
    mcell.z = mcell.z * cosd(cell_beta - 90);

elseif moltype==260 %r26001_Cytosine_Barker.xyz
    cell_a=13.041;
    cell_b=9.494;
    cell_c=3.815;
    cell_beta=90;

    DA1x=0; %screw diad axis [0 y 1/2]
    DA1z=1/2;

    DA2y=1/4; %screw diad axis [x 1/4 1/2]
    DA2z=1/2;

    DA3x=1/2; %screw diad axis [1/2 y 1/2]
    DA3z=1/2;

    mc1 = mc0; 
    mc1.x = -(mc1.x-DA1x) + DA1x;
    mc1.y = mc1.y + 1/2;
    mc1.z = -(mc1.z-DA1z) + DA1z;

    mc2 = mc0; 
    mc2.x = mc2.x + 1/2;
    mc2.y = -(mc2.y-DA2y) + DA2y;
    mc2.z = -(mc2.z-DA2z) + DA2z;

    mc3 = mc2; 
    mc3.x = -(mc3.x-DA3x) + DA3x;
    mc3.y = mc3.y + 1/2;         
    mc3.z = -(mc3.z-DA3z) + DA3z;

    mcell=mc0;
    mcell=mergemol(mcell,mc1);
    mcell=mergemol(mcell,mc2);
    mcell=mergemol(mcell,mc3);

    mcell.x=mcell.x*cell_a;  %denormalizing coordinates
    mcell.y=mcell.y*cell_b;
    mcell.z=mcell.z*cell_c;

    mcell.x = mcell.x + mcell.z.*cosd(cell_beta);
    mcell.z = mcell.z * cosd(cell_beta - 90);

elseif moltype==261 %
    cell_a=7.801;
    cell_b=9.844;
    cell_c=7.683;
    cell_beta=99+42/60;

    DA1x=0; %screw diad axis [0 y 1/2]
    DA1z=1/2;

    DA2y=1/4; %screw diad axis [x 1/4 1/2]
    DA2z=1/2;

    DA3x=1/2; %screw diad axis [1/2 y 1/2]
    DA3z=1/2;

    mc1 = mc0; 
    mc1.x = -(mc1.x-DA1x) + DA1x;
    mc1.y = mc1.y + 1/2;
    mc1.z = -(mc1.z-DA1z) + DA1z;

    mc2 = mc0; 
    mc2.x = mc2.x + 1/2;
    mc2.y = -(mc2.y-DA2y) + DA2y;
    mc2.z = -(mc2.z-DA2z) + DA2z;

    mc3 = mc2; 
    mc3.x = -(mc3.x-DA3x) + DA3x;
    mc3.y = mc3.y + 1/2;         
    mc3.z = -(mc3.z-DA3z) + DA3z;

    mcell=mc0;
    mcell=mergemol(mcell,mc1);
    mcell=mergemol(mcell,mc2);
    mcell=mergemol(mcell,mc3);

    mcell.x=mcell.x*cell_a;  %denormalizing coordinates
    mcell.y=mcell.y*cell_b;
    mcell.z=mcell.z*cell_c;

    mcell.x = mcell.x + mcell.z.*cosd(cell_beta);
    mcell.z = mcell.z * cosd(cell_beta - 90);
elseif moltype==12 %dUrd_Rahman
    cell_a=7.91;
    cell_b=6.710;
    cell_c=18.77;
    cell_beta=96.6;



    mcell=mc0;


    mcell.x=mcell.x*cell_a;  %denormalizing coordinates
    mcell.y=mcell.y*cell_b;
    mcell.z=mcell.z*cell_c;

    mcell.x = mcell.x + mcell.z.*cosd(cell_beta);
    mcell.z = mcell.z * cosd(cell_beta - 90);
elseif moltype==140 %rAdo_Klooster
    cell_a=4.78858;%Ang
    cell_b=10.2402;%Ang
    cell_c=11.7722;%Ang
    cell_beta=99.592;%deg

    mcell=mc0;

    mcell.x=mcell.x*cell_a;  %denormalizing coordinates
    mcell.y=mcell.y*cell_b;
    mcell.z=mcell.z*cell_c;

    mcell.x = mcell.x + mcell.z.*cosd(cell_beta);
    mcell.z = mcell.z * cosd(cell_beta - 90);
elseif moltype==9 %dCyd_Subraminain
    cell_a=6.561;%Ang
    cell_b=17.659;%Ang
    cell_c=5.125;%Ang
    cell_beta=108.08;%deg

    mcell=mc0;

    mcell.x=mcell.x*cell_a;  %denormalizing coordinates
    mcell.y=mcell.y*cell_b;
    mcell.z=mcell.z*cell_c;

    mcell.x = mcell.x + mcell.z.*cosd(cell_beta);
    mcell.z = mcell.z * cosd(cell_beta - 90);
elseif moltype==16 %dGuo_Haschmeyer
    %P22_1 2_1
    cell_a=5.14;%Ang
    cell_b=19.11;%Ang
    cell_c=23.66;%Ang
    cell_beta=90;%deg

    mcell=mc0;

    mcell.x=mcell.x*cell_a;  %denormalizing coordinates
    mcell.y=mcell.y*cell_b;
    mcell.z=mcell.z*cell_c;

    mcell.x = mcell.x + mcell.z.*cosd(cell_beta);
    mcell.z = mcell.z * cosd(cell_beta - 90);
end


mcell = createbondtable(mcell);

mcell.desc=[mc0.desc '.out'];
if fl_writefile
    savemol(workdir,mcell);
end

if fl_repcells
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

    mcellf=mcell;
    mcellf=mergemol(mcellf,mcell3);
    if 0
        mcellf=mergemol(mcellf,mcell2);
        mcellf=mergemol(mcellf,mcell5);
    end
    if 1
        mcellf=mergemol(mcellf,mcell4);
        mcellf=mergemol(mcellf,mcell6);
%        mcellf=mergemol(mcellf,mcell7);
%        mcellf=mergemol(mcellf,mcell8);
    end

    mcellf = createbondtable(mcellf);
    plotmol(mcellf,'k',0,0);

    mcellf.desc=[mc0.desc '.outf'];
    if fl_writefile
        savemol(workdir,mcellf);
    end

else
    plotmol(mcell);
end



plot3([0 cell_a cell_a 0],[-cell_b/2 -cell_b/2 -cell_b/2 -cell_b/2],[0 0 cell_c cell_c],'k');
plot3([0 cell_a cell_a 0],[-cell_b/4 -cell_b/4 -cell_b/4 -cell_b/4],[0 0 cell_c cell_c],'k');
plot3([0 cell_a cell_a 0],[ cell_b/4  cell_b/4  cell_b/4  cell_b/4],[0 0 cell_c cell_c],'k');
plot3([0 cell_a cell_a 0],[ cell_b/2  cell_b/2  cell_b/2  cell_b/2],[0 0 cell_c cell_c],'k');

plot3([0 0],[-cell_b/2  cell_b/2],[0 0],'k');
plot3([cell_a cell_a],[-cell_b/2  cell_b/2],[0 0],'k');
plot3([0 0],[-cell_b/2  cell_b/2],[cell_c cell_c],'k');
plot3([cell_a cell_a],[-cell_b/2  cell_b/2],[cell_c cell_c],'k');



