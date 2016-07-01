%add guanosine base to selected sugar conformations and make deoxyguanine anti and syn conformations
%
% Version 1.0    
% Last modified  R O Zhurakivsky 2008-05-19
% Created        R O Zhurakivsky 2006-?-?

format compact
clear
global pind
global flPlot

atomsind
flPlot=1;

workdbname=[CD.dbdir filesep 'r8_g_or.mat']
workname='rg01';
odir=[CD.xyzdir filesep workname];
moltype=16;
gtemplname=[workname '_templ.gjf']
fullgtemplname=[CD.templatesdir filesep gtemplname]


if 0
workdbname1='r11.mat'
load(workdbname1,'workdb')
set1desc=cell(0);
recnum=numel(workdb);
for i=1:recnum
%  if workdb(i).new=='Y'
%    set1desc(end+1)={workdb(i).filename}; %names of result files I have
    set1desc(end+1)={workdb(i).prop.sdesc}; %names of result files I have
%  end
end
end
set1desc={};

afile = 'guanosine.xyz'
fullafilename=[CD.templatesdir filesep afile]
aH9 = 1; %indexes of glycosidic linkage atoms %H9
aN9 = 2; %N9
aC4 = 16;
aN7 = 5;

dist_glyc = (1.4595+1.4682)/2; %the mean value between AaaaA & AaaaS dAde conformations (rf01)

mc0=loadmolxyz(fullafilename);

load(workdbname,'workdb')

recnum=numel(workdb);
for i=1:recnum
for j=1:2

    ms0=workdb(i);

    if j==1
        mcyd_desc = [ms0.desc 'S'];
        changefi = 0;
    else
        mcyd_desc = [ms0.desc 'A'];
        changefi = 1;
    end

    [iprev,inext,sdesc_prev,sdesc_next] = checksimdesc(mcyd_desc,set1desc);
    if ~isempty(strcmpcellar(set1desc,mcyd_desc)) | iprev | inext
       continue
    end


    ['Loaded conformation: ' ms0.prop.sdesc]
%plotmol(ms0);

    maxind=0;
%    maxind=ms0.atomnum;
    mc0.atomnum = length(mc0.labels);
    mc0.ind = ((maxind+1):(maxind+mc0.atomnum))';
    maxind=maxind+mc0.atomnum;


    iC1  = ms0.ind(find(find(strcmp(pind.labels,'pC1')) ==ms0.pind)); 
    iH11 = ms0.ind(find(find(strcmp(pind.labels,'pH11'))==ms0.pind)); 
    iH12 = ms0.ind(find(find(strcmp(pind.labels,'pH12'))==ms0.pind)); 


    %vector C1H11 of sugar
    vect_sug = createvect(ms0,iC1,ms0,iH11);
    cos_vect_sug = vect_sug/norm(vect_sug);

    %movement vector of N1 atom of cytidine
    vect_mov_cyt = createvect(mc0,aN9,ms0,iC1)+cos_vect_sug*dist_glyc;

%    plotmol(mc0,'r');

    mc1 = mc0;
    mc1.x = mc0.x + vect_mov_cyt(1);
    mc1.y = mc0.y + vect_mov_cyt(2);
    mc1.z = mc0.z + vect_mov_cyt(3);
%    plotmol(mc1,'k');

    %vector H1-N1 of cytidine
    vect_cyt = createvect(mc1,aN9,mc1,aH9);

    %rotates cytidine so that N1H1 to be along C1'H11'
    mc2 = rotvect3(-vect_sug, vect_cyt, mc1, aN9);
    %rotates cytidine to create syn- conformation
%    plotmol(mc2,'g');
    fi = p2pangle(createplane(ms0,iH12,ms0,iC1,ms0,iH11),createplane(mc2,aN7,mc2,aN9,mc2,aC4));
    %determine is H12' above or under C2N1C6 plane and correct rotation angle according
    v1=createvect(mc2,aN9,mc2,aN7);
    v2=createvect(mc2,aN9,mc2,aC4);
    v3=createvect(ms0,iC1,ms0,iH12);
    r = dot(cross(v1,v2),v3);
    if r<0
        fi=-fi;
    end
    if changefi
        fi=fi+pi;
    end


%    vect_cyt = createvect(mc2,aN1,mc2,aH1);
    mc3 = rotline3(vect_sug,mc2,fi,aN9);

    ms1 = delatom( ms0, iH11 );
    mc3 = delatom( mc3, aH9 );
    drawbond(ms1,iC1,mc3,aN9);

%    plotmol(mc2,'k');
%    plotmol(ms0);

    molmol=mergemol(ms1,mc3);
    molmol.desc = strcat(workname,'_',mcyd_desc);
    molmol = createbondtable(molmol);
    molmol = identmol(molmol,moltype); 

    if exist(odir)~=7
        mkdir(odir);
    end
    savemol(odir,molmol,0);
%    [xxx,order]=sortrows(molmol.pind);
%    savemolgs(odir,molmol,4,order,fullgtemplname); %Gaussian with ZMT
order=[];
    savemolgs(odir,molmol,3,order,fullgtemplname); %Gaussian with ZMT

plotmol(molmol);

%pause
clf

end
end

