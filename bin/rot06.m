%tests
%ZMT to XYZ conversion
%
% Version 1.0    
% Last modified  R O Zhurakivsky 2007-08-04
% Created        R O Zhurakivsky 2005-09-?

atomsind
indir=[CD.xyzdir filesep 'testzmt'];

clear mz0
format compact
grid on

odir='out4';

sfiles = dir(strcat(indir,filesep,'*.zmt'));

numfiles = size(sfiles,1);

maxind=0;
blen=zeros(numfiles,1);
for i=1:numfiles
    delete(gca);

    dlm=strfind(sfiles(i).name,'.');
    fname = sfiles(i).name(1:(dlm-1));

    ffname = strcat(indir,filesep,sfiles(i).name); 

%    [mz0.labels,mz0.iR,mz0.R,mz0.ialfa,mz0.alfa,mz0.ibeta,mz0.beta,]=textread(ffname,'%c%d%f%d%f%d%f');
    [mz0.labels,mz0.iR,mz0.R,mz0.ialfa,mz0.alfa,mz0.ibeta,mz0.beta]=textread(ffname,'%c%d%f%d%f%d%f')';
    mz0.iR=uint16(mz0.iR);
    mz0.ialfa=uint16(mz0.ialfa);
    mz0.ibeta=uint16(mz0.ibeta);

    mz0.atomnum = uint16(length(mz0.labels));
    mz0.ind = ((maxind+1):(maxind+mz0.atomnum))';
    mz0.desc = fname;

    mz0=zmt2xyz(mz0,1:mz0.atomnum);
        
    strcat('Loaded file: ',ffname)

    plotmol(mz0);
    
    mz0.pind=1:mz0.atomnum;
    savemol(odir,mz0,0); 

    
%    'Press any key.'
%    pause
end

