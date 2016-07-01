%tests
%make all possible conformations by rotating around ordinary bonds and
%create NWChem input file
%
% Version 1.0    
% Last modified  R O Zhurakivsky 2007-03-17
% Created        R O Zhurakivsky 2005-09-?

atomsind
global pind
indir=[CD.xyzdir filesep '01'];

clear ms0
format compact
grid on

odir=[CD.xyzdir filesep 'out5'];

sfiles = dir(strcat(indir,filesep,'*.xyz'));

numfiles = size(sfiles,1);

maxind=0;
blen=zeros(numfiles,1);
for i=1:numfiles
    delete(gca);

    dlm=strfind(sfiles(i).name,'.');
    fname = sfiles(i).name(1:(dlm-1));

    ffname = strcat(indir,filesep,sfiles(i).name); 
	ms0=loadmolxyz(ffname,fname,maxind);

    ms0 = createbondtable(ms0);
    ms0 = identmol(ms0);

    strcat('Loaded file: ',ffname)

%    [y,order]=sortrows(ms0.pind); 
    ms0 = createzmt(ms0);

    iO2 = ms0.ind(find(find(strcmp(pind.labels,'pO2'))==ms0.pind)); 
    iC2 = ms0.ind(find(find(strcmp(pind.labels,'pC2'))==ms0.pind));

    ind = findallrotates(ms0,iO2,iC2);
    

    ms0=zmt2xyz(ms0);
    plotmol(ms0);
    savemolnw(odir,ms0,3); 

    ms1 = rotaroundbond(ms0,ind,120);
    ms1.desc = strcat(ms1.desc,'-0001');
    plotmol(ms1);
    savemolnw(odir,ms1,3); 

    ms2 = rotaroundbond(ms0,ind,240);
    ms2.desc = strcat(ms2.desc,'-0002');
    plotmol(ms2);
    savemolnw(odir,ms2,3); 


return    
%    'Press any key.'
%    pause
end

