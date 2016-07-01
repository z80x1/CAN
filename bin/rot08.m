%tests
%calculate molecule properties
%
% Version 1.0    
% Last modified  R O Zhurakivsky 2007-03-17
% Created        R O Zhurakivsky 2005-09-?

clear 
format compact
%grid on

atomsind
global pind
%indir=CD.xyzdir filesep '01'];
indir='..\work\_yevgen\060530\'


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
plotmol(ms0); 
%    ms0 = identmol(ms0);

    strcat('Loaded file: ',ffname)

%    ms0 = createzmt(ms0);
%    ms0=zmt2xyz(ms0);

%    ms0 = calcproperties(ms0);

    
%    dbrib(i)=ms0;

disp('Press any key.')
    pause
end
%save 'dbrib.mat' dbrib
