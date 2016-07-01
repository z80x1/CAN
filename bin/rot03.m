%calculates H-bonds distance map   
%
% Version 1.0    
% Last modified  R O Zhurakivsky 2007-03-17
% Created        R O Zhurakivsky 2005-09-?

global flPlot
%flPlot=1;

%load atoms indexes
clear

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

dist_glik = 1.4743; %the mean value of 25 optimized cytidine conformations

minHlen = 1.1;
maxHlen = 2.5;

indir=[CD.xyzdir filesep 'drib'];
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


%files = ls(strcat(indir,filesep,'*.xyz')); %M7
sfiles = dir(strcat(indir,filesep,'*.xyz'));

numfiles = size(sfiles,1);
distmap=zeros(numfiles,360/5,17);
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
    vect_mov_cyt = createvect(ms0,aC1_,mc0,aN1)+cos_vect_sug*dist_glik;


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
%    [fi,r]
    if r>0
        fi=-fi;
    end
    mc2 = rotline3(vect_sug,mc2,fi,aN1);


    ms0 = delatom( ms0, aH11_ );
    mc2 = delatom( mc2, aH1 );
    drawbond(ms0,aC1_,mc2,aN1);

    plotmol(ms0);

        distmapstr='';
%        distmap=zeros(360/5,17);
        fi_min = 0;
        fi_step = 5*pi/180;
        fi_max = 2*pi;
        j=0;
        for fi=fi_min:fi_step:fi_max,
            if fi~=0
                mc3 = rotline3(vect_sug,mc2,fi,aN1);
            else
                mc3=mc2;
            end

            dH5 = adist([ms0 mc3],aH5_,aO2);
            dO2O5 = adist([ms0 mc3],aO5_,aO2);
            angH5 = valang([ms0 ms0 mc3],aO5_,aH5_,aO2);

            dH3 = adist([ms0 mc3],aH3_,aO2);
            dO2C3 = adist([ms0 mc3],aC3_,aO2);
            angH3 = valang([ms0 ms0 mc3],aC3_,aH3_,aO2);

            dH2 = adist([ms0 mc3],aH2_,aO2);
            dO2C2 = adist([ms0 mc3],aC2_,aO2);
            angH2 = valang([ms0 ms0 mc3],aC2_,aH2_,aO2);

            dH1 = adist([ms0 mc3],aH12_,aO2);
            dO2C1 = adist([ms0 mc3],aC1_,aO2);
            angH1 = valang([ms0 ms0 mc3],aC1_,aH12_,aO2);


            dH6H5 = adist([ms0 mc3],aH5_,aH6);
            dH6H3 = adist([ms0 mc3],aH3_,aH6);
            dH6H2 = adist([ms0 mc3],aH2_,aH6);
            dH6H1 = adist([ms0 mc3],aH12_,aH6);


            nbond = 0;

            if (dO2O5 < 3) | (dO2C3 < 3) | (dO2C1 < 3) | (dO2C1 < 3)
              j=j+1;            
              distmap(i,j,:) = [fi*180/pi, dH5,dO2O5,angH5, dH3,dO2C3,angH3, dH2,dO2C2,angH2, dH1,dO2C1,angH1,  dH6H5,dH6H3,dH6H2,dH6H1];
              ss=sprintf('%3.0f   %4.2f %4.2f %3.0f    %4.2f %4.2f %3.0f    %4.2f %4.2f %3.0f    %4.2f %4.2f %3.0f      %4.2f   %4.2f   %4.2f   %4.2f',...
fi*180/pi, dH5,dO2O5,angH5, dH3,dO2C3,angH3, dH2,dO2C2,angH2, dH1,dO2C1,angH1,  dH6H5,dH6H3,dH6H2,dH6H1);
              distmapstr=strvcat(distmapstr,ss);

%              m1=mergemol(ms0,mc3);
%              m1.desc = strcat(f1,'_',fname,'_',num2str(round(fi*180/pi),'%03d'));
%              savemol(outdir,m1);
            end
        end
%        ['  fi          dH5''O2    dO5''O2   angH5       dH3''O2    dC3''O2   angH3      dH2''O2    dC2''O2   angH2''      dH1''O2    dC1''O2   angH1']
        ['fi dH5''O2 dO5''O2 angH5  dH3''O2 dC3''O2 angH3  dH2''O2 dC2''O2 angH2''  dH1''O2 dC1''O2 angH1    dH6H5 dH6H3 dH6H2 dH6H1']
        distmapstr
%        indmin =  find(distmap(i,1:j,9)==min(distmap(i,1:j,9))); %dC2'O2 min
%        distmapstr(indmin,:)
        if nbond==0
            plotmol(mc2,'r');
        end            
%    'Press any key.'
    pause
end
