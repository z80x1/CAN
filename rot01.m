%   finds sugar conformations that have only O5'H...O2 H-bond and haven't 
%O5'H...O4' one
%
% Version 1.0    
% Last modified  R O Zhurakivsky 2007-03-17
% Created        R O Zhurakivsky 2005-09-?

global flPlot
%flPlot=0;

clear
%load atoms indexes
%atoms indexes
aH11_ = 23;
aH12_ = 14;
aC1_  = 15;
aC2_ = 16;
aC3_ = 17;
aH2_ = 25;
aH3_ = 26;
aO3_ = 22;
aO4_ = 21;
aO5_ = 20;
aH5_ = 30;

aN1 = 2;
aH1 = 1;
aO2 = 13;
aC2 = 12;
aC6 = 3;
aH6 = 4;

atomsind

hbond_ind1 = [aO5_ aH5_ aO4_];
hbond_ind2 = [aO5_ aH5_ aO2];
%hbond_ind2 = [aO5_ aH5_ aO4_];
%hbond_ind2 = [aO5_ aH5_ aO3_];
%hbond_ind2 = [aO5_ aH3_ aO3_];


dist_glyc = 1.4743; %the mean value of 25 optimized cytidine conformations

indir=[CD.xyzdir filesep  'drib'];
file1='cyt_base.xyz';
outdir=[CD.xyzdir filesep 'out'];

format compact
grid on
maxind=0;

mc0=loadmolxyz(file1,'',maxind);

maxind=maxind+mc0.atomnum;

%dlm=strfind(file1,'.'); %M7
dlm=strfind(file1,'.');
f1 = file1(1:(dlm-1));


%files = ls(strcat(indir,'/*.xyz')); %M7
sfiles = dir(strcat(indir,filesep,'*.xyz'));

numfiles = size(sfiles,1);
for i=1:numfiles

    delete(gca);

%    dlm=strfind(files(i,:),'.'); %M7
%    fname = files(i,1:(dlm-1));  %M7
    dlm=strfind(sfiles(i).name,'.');
    fname = sfiles(i).name(1:(dlm-1));

    ffname = strcat(indir,filesep,sfiles(i).name); 
	ms0=loadmolxyz(ffname,'',maxind);

    strcat('Loaded file: ',ffname)

    %vector C1H11 of sugar
    vect_sug = createvect(ms0,aH11_,ms0,aC1_);
    cos_vect_sug = vect_sug/norm(vect_sug);

    %movement vector of N1 atom of cytidine
    vect_mov_cyt = createvect(ms0,aC1_,mc0,aN1)+cos_vect_sug*dist_glyc;


    mc1 = mc0;
    mc1.x = mc0.x + vect_mov_cyt(1);
    mc1.y = mc0.y + vect_mov_cyt(2);
    mc1.z = mc0.z + vect_mov_cyt(3);
    %plotmol(mc1,'g');


    %vector H1-N1 of cytidine
    vect_cyt = createvect(mc1,aN1,mc1,aH1);

    %rotates cytidine so that N1H1 to be along C1'H11'
    mc2 = rotvect3(vect_sug, vect_cyt, mc1, aN1);
    %rotates cytidine to create syn- conformation
    %plotmol(mc2,'g');
    fi = p2pangle(createplane(ms0,aH12_,ms0,aC1_,ms0,aH11_),createplane(mc2,aC6,mc2,aN1,mc2,aC2));
    %determine is H12' above or under C2N1C6 plane and correct rotation angle according
    v1=createvect(mc2,aN1,mc2,aC6);
    v2=createvect(mc2,aN1,mc2,aC2);
    v3=createvect(ms0,aC1_,ms0,aH12_);
    r = dot(cross(v1,v2),v3);
    if r>0
        fi=-fi;
    end
    mc2 = rotline3(vect_sug,mc2,fi,aN1);


    ms0 = delatom( ms0, aH11_ );
    mc2 = delatom( mc2, aH1 );
    drawbond(ms0,aC1_,mc2,aN1);

    plotmol(ms0);

    d0 = adist(ms0,hbond_ind1(1),hbond_ind1(3));
    ang0 = valang(ms0,hbond_ind1(1),hbond_ind1(2),hbond_ind1(3));
    [d0 ang0] %#ok
    if (d0 <=3 ) && (ang0 >= 90)
        strcat('H-bond O5H5...O4 exists.')
        plotmol(mc2,'r');

        if 0
            m1=mergemol(ms0,mc2);
            m1.desc = strcat(f1,'_',fname);
            savemol(outdir,m1);
        end
    else

        fi_min = 5*pi/180;
        fi_step = 5*pi/180;
        fi_max = 2*pi;
        nbond = 0;
        for fi=fi_min:fi_step:fi_max,
            mc3 = rotline3(vect_sug,mc2,fi,aN1);
%            mc3=mc2;
            m1=mergemol(ms0,mc3);
            d = adist([ms0 mc3],hbond_ind2(1),hbond_ind2(3));
            ang = valang([ms0 ms0 mc3],hbond_ind2(1),hbond_ind2(2),hbond_ind2(3));
%            d = adist(ms0,hbond_ind2(1),ms0,hbond_ind2(3));
%            ang = valang(ms0,hbond_ind2(1),ms0,hbond_ind2(2),ms0,hbond_ind2(3));
            if (d<=3) && (ang >= 90)
                nbond = nbond +1;
                [fi*180/pi d ang]  %#ok
                plotmol(mc3,'r');
%                m1.desc = strcat(f1,'_',fname,'_',num2str(round(fi*180/pi),'%03d'));
%                savemol(outdir,m1);
            end
        end
        if nbond==0
            plotmol(mc2,'r');
        end            
    end
%    'Press any key.'
   pause
end
