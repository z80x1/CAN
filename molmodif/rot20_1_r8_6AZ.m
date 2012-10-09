%add azapirimidine ring to selected sugar conformations and make anti and syn conformations
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
workname='rb04';
odir=[CD.xyzdir filesepworkname];
moltype=11;
gtemplname=[workname '_templ.gjf']
fullgtemplname=[CD.templatesdir filesep gtemplname]


workdbname1='r11.mat' %database with already obtained conformations
load(workdbname1,'workdb')
set1desc=cell(0);
recnum=numel(workdb);
for i=1:recnum
%  if workdb(i).new=='Y'
%    set1desc(end+1)={workdb(i).filename}; %names of result files I have
    set1desc(end+1)={workdb(i).prop.sdesc}; %names of result files I have
%  end
end


afile = 'azacyt_base.xyz'
aH1 = 9; %indexes of glycosidic linkage atoms %H1
aN1 = 2; %N1
aC2 = 7;
aN6 = 3;

dist_glyc = 1.476; %the mean value of 82 optimized AC6 conformations (r11.050813)


mc0=loadmolxyz(afile);

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
    vect_mov_cyt = createvect(mc0,aN1,ms0,iC1)+cos_vect_sug*dist_glyc;

%    plotmol(mc0,'r');

    mc1 = mc0;
    mc1.x = mc0.x + vect_mov_cyt(1);
    mc1.y = mc0.y + vect_mov_cyt(2);
    mc1.z = mc0.z + vect_mov_cyt(3);
%    plotmol(mc1,'k');

    %vector H1-N1 of cytidine
    vect_cyt = createvect(mc1,aN1,mc1,aH1);

    %rotates cytidine so that N1H1 to be along C1'H11'
    mc2 = rotvect3(-vect_sug, vect_cyt, mc1, aN1);
    %rotates cytidine to create syn- conformation
%    plotmol(mc2,'g');
    fi = p2pangle(createplane(ms0,iH12,ms0,iC1,ms0,iH11),createplane(mc2,aN6,mc2,aN1,mc2,aC2));
    %determine is H12' above or under C2N1C6 plane and correct rotation angle according
    v1=createvect(mc2,aN1,mc2,aN6);
    v2=createvect(mc2,aN1,mc2,aC2);
    v3=createvect(ms0,iC1,ms0,iH12);
    r = dot(cross(v1,v2),v3);
    if r<0
        fi=-fi;
    end
    if changefi
        fi=fi+pi;
    end


%    vect_cyt = createvect(mc2,aN1,mc2,aH1);
    mc3 = rotline3(vect_sug,mc2,fi,aN1);

    ms1 = delatom( ms0, iH11 );
    mc3 = delatom( mc3, aH1 );
    drawbond(ms1,iC1,mc3,aN1);

%    plotmol(mc2,'k');
%    plotmol(ms0);

    mcyd=mergemol(ms1,mc3);
    mcyd.desc = strcat(workname,'_',mcyd_desc);
    mcyd = createbondtable(mcyd);
    mcyd = identmol(mcyd,moltype); 

    if exist(odir)~=7
        mkdir(odir);
    end
    savemol(odir,mcyd,0);
%not needed
%    [xxx,order]=sortrows(mcyd.pind);
order=[];
    savemolgs(odir,mcyd,4,order,fullgtemplname); %Gaussian with ZMT

plotmol(mcyd);

%pause
clf

end
end

