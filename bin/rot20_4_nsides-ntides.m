%add P contaning residue to selected sugar conformations and make
%nucleotide sugar residue

%merge molA with molB so that iA_2 and iB_2 atoms to be at dist_link
%distance. iA_1 & iB_1 atoms will be deleted.
%
% Version 1.0    
% Last modified  R O Zhurakivsky 2008-05-19
% Created        R O Zhurakivsky 2006-?-?


%molA - first molecule
%molB - second molecule
%molZ - resulting molecule

format compact
clear
global pind
global flPlot

atomsind
flPlot=1;

moltype0=12  %#ok
theory='dftV2'  %#ok
onlyoriginal=1;  % process db with only original conformations


workname='r53003';
odir=[CD.xyzdir filesep workname]; %output directory
moltype=530; %type of resulting molecule
gtemplname=[workname '_templ.gjf']
fullgtemplname=[CD.templatesdir filesep gtemplname]
sortmode=2  %(0/1/2 - without sorting/sort by sdesc/energy)

set1desc={};

molAfname = 'Pres.xyz' %mol_B
fullmolAfname=[CD.templatesdir filesep molAfname]

if ~strcmp(theory,'dft')
  theorystr = ['_' theory];
else
  theorystr = '';
end
workdbname=[CD.dbdir filesep 'r' int2str(moltype0) '_g' theorystr];
if onlyoriginal
    templ='_or';
    workdbname = [workdbname templ];
end
workdbname=[workdbname '.mat']  %#ok  %mol_A  DB

iA_1=2; %indexes of OH group of P residue needed to be connected to sugar (A#1) (aHP4)
iA_2=3; %(A#2) (aOP4)
iA_3=5; %(A#3) (aOP3 -  O atom without hydrogen) 
iA_2delete = [1]; %indexes of additional atoms to be deleted

%aH9 = 1; %indexes of glycosidic linkage atoms %H9
%aN9 = 2; %N9
%aC4 = 16;
%aN7 = 5;

dist_link = 1.63; %wanted distance between P residue and sugar
dist_tor = 347.8756/180*pi; %wanted torsion value of iA_3-iA_2-iB_2-iB_3


molA0=loadmolxyz(fullmolAfname);

if numel(iA_2delete) > 1
   error('deleting of multiply atoms is not supported');
end
for i=iA_2delete
    molA0 = delatom( molA0, i );
    %after atoms deleting atom numbering is changed :(
    iA_1=iA_1-(iA_1>i);
    iA_2=iA_2-(iA_2>i);
    iA_3=iA_3-(iA_3>i);
end
    
load(workdbname,'workdb')

recnum=numel(workdb);
for i=1:recnum

    if isfield(workdb(i),'GO') & isfield(workdb(i).GO,'energy')
        GOenergy(i) = workdb(i).GO.energy;
    else
        sortbyenergy=0;
    end
    sdesc(i) = {workdb(i).prop.sdesc};
end

if sortmode==1 %sort by sdesc
  [xxx,ind]=sort(sdesc);
elseif sortmode==2 %sort by energy
  [xxx,ind]=sort(GOenergy);
else %without sorting
  ind=1:recnum;
end

fileind=0;
for i=ind
%for j=1:2  % second cycle is for creating syn/anti conformations

    molB0=workdb(i);
    molZ_desc = molB0.desc;

    [iprev,inext,sdesc_prev,sdesc_next] = checksimdesc(molZ_desc,set1desc);
    if ~isempty(strcmpcellar(set1desc,molZ_desc)) | iprev | inext
       continue
    end

    ['Loaded conformation: ' molB0.prop.sdesc]
    fileind=fileind+1;
%plotmol(molB0);

if 0
    iB_1 = molB0.ind(find(find(strcmp(pind.labels,'pH53'))==molB0.pind)); %atom that will be deleted
    iB_2  = molB0.ind(find(find(strcmp(pind.labels,'pO5')) ==molB0.pind)); %inner atom of future links
    iB_3  = molB0.ind(find(find(strcmp(pind.labels,'pC5')) ==molB0.pind)); %third atom in chain
else
%    molB0=molZ;
    iB_1 = molB0.ind(find(find(strcmp(pind.labels,'pH32'))==molB0.pind)); %atom that will be deleted
    iB_2  = molB0.ind(find(find(strcmp(pind.labels,'pO3')) ==molB0.pind)); %inner atom of future links
    iB_3  = molB0.ind(find(find(strcmp(pind.labels,'pC3')) ==molB0.pind)); %third atom in chain
end

    maxind=0;
%    maxind=molB0.atomnum;
    molA0.atomnum = length(molA0.labels);
    molA0.ind = ((maxind+1):(maxind+molA0.atomnum))';
    maxind=maxind+molA0.atomnum;



    %vector #2#1 of second molecule
    vect_B = createvect(molB0,iB_2,molB0,iB_1);
    cos_vect_B = vect_B/norm(vect_B);

    %movement vector of #2 atom of first molecule
    vect_mov_A = createvect(molA0,iA_2,molB0,iB_2)+cos_vect_B*dist_link;

%    plotmol(molA0,'r');

    molA1 = molA0;
    molA1.x = molA0.x + vect_mov_A(1);
    molA1.y = molA0.y + vect_mov_A(2);
    molA1.z = molA0.z + vect_mov_A(3);
%    plotmol(molA1,'k');

    %vector #2-#1 of first molecule
    vect_A = createvect(molA1,iA_2,molA1,iA_1);

    %rotates first molecule so that #2#1 vector to be along #2#1 vector of second molecule
    molA2 = rotvect3(-vect_B, vect_A, molA1, iA_2);

    %constrain iA_3-iA_2-iB_2-iB_3 torsion
    %%rotates cytidine to create syn- conformation
%    plotmol(molA2,'g');
    fi = p2pangle(createplane(molB0,iB_1,molB0,iB_2,molB0,iB_3),createplane(molA2,iA_1,molA2,iA_2,molA2,iA_3))+pi		
    %determine is H12' above or under C2N1C6 plane and correct rotation angle according
%    v1=createvect(molA2,iA_2,molA2,iA_1);
%    v2=createvect(molA2,iA_2,molA2,iA_3);
%    v3=createvect(molB0,iB_2,molB0,iB_3);
%    r = dot(cross(v1,v2),v3);
%    if r<0
%        fi=-fi
%    end
%    if changefi
%       fi=fi+pi;
%    end

%    vect_A = createvect(molA2,aN1,molA2,aH1);
    molA3 = rotline3(vect_B,molA2,dist_tor-fi,iA_2);

    molB1 = delatom( molB0, iB_1 );
    molA3 = delatom( molA3, iA_1 );
    drawbond(molB1,iB_2,molA3,iA_2);

%    plotmol(molA2,'k');
%    plotmol(molB0);

    molZ=mergemol(molB1,molA3);
%    molZ.desc = strcat(workname,'_',sprintf('%03d',fileind),'_',molZ_desc);
    molZ.desc=[workname,'#r',int2str(moltype),'_',sprintf('%03d',fileind),'_',molB0.prop.sdesc];

    molZ = createbondtable(molZ);
%    molZ = identmol(molZ,moltype); 

    if exist(odir)~=7
        mkdir(odir);
    end
    savemol(odir,molZ,0);

order=[];
    savemolgs(odir,molZ,3,order,fullgtemplname); %Gaussian with ZMT
%plotmol(molZ,'b',0);

%pause
%clf

%end
end

