%   calculates average C1'N1 bond length over optimized cytidine conformations
% in 'cyt' directory
%
% Version 1.0    
% Last modified  R O Zhurakivsky 2007-03-17
% Created        R O Zhurakivsky 2005-09-?


indir=[CD.xyzdir filesep 'cyt'];

clear ms0
format compact
grid on


blen=zeros(numfiles,1);

sfiles = dir(strcat(indir,filesep,'*.xyz'));

numfiles = size(sfiles,1);
for i=1:numfiles

    delete(gca);

    dlm=strfind(sfiles(i).name,'.');
    fname = sfiles(i).name(1:(dlm-1));

    ffname = strcat(indir,filesep,sfiles(i).name); 
	ms0=loadmolxyz(ffname,'',maxind);

    strcat('Loaded file: ',ffname)


%    plotmol(ms0);
    blen(i)=adist(ms0,14,27);

%    'Press any key.'
%    pause
end

blen
mean(blen)