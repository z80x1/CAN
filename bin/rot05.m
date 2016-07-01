%tests
%Z-Matrix creation
%ordered by molident
%
% Version 1.0    
% Last modified  R O Zhurakivsky 2007-03-17
% Created        R O Zhurakivsky 2005-09-?

moltype=8 %1'2'-deoxyribose

atomsind

indir=[CD.xyzdir filesep 'dryb.xyz'];

clear ms0
format compact
grid on

odir=[CD.xyzdir filesep 'out3'];

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
    ms0 = identmol(ms0,moltype); 

    strcat('Loaded file: ',ffname)

%    plotmol(ms0);
    
%    [y,i]=sortrows(ms0.pind);

%    ms0 = createzmt(ms0,i');

    savemol(odir,ms0,2); %ZMT
%    'Press any key.'
%    pause
end

