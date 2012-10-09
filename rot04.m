%tests
%%Z-Matrix creation
%
% Version 1.0    
% Last modified  R O Zhurakivsky 2007-03-17
% Created        R O Zhurakivsky 2005-09-?

atomsind
indir=[CD.xyzdir filesep 'cyt'];

clear ms0
format compact
grid on

odir=[CD.xyzdir filesep 'cyt.out'];

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

    strcat('Loaded file: ',ffname)


    plotmol(ms0);
    

    savemol(odir,ms0,1);

%    'Press any key.'
    pause
end

