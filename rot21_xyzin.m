%rot21: read xyz files, identify them and output xyz or zmt files in internal coordinates
%
% Version 1.0    
% Last modified  R O Zhurakivsky 2007-03-17
% Created        R O Zhurakivsky 2006-?-?

moltype=12
%indir=[CD.xyzdir '/' 'rc17/eaacs']
indir=[CD.datadir filesep 'g.out/rc/rc17/eaacs']
odir=[indir filesep 'out_geom']


global pind

clear ms0 workdb desc filenames
format compact
atomsind


sfiles = dir(strcat(indir,filesep,'*.xyz'));
numfiles = size(sfiles,1);

for i=1:numfiles

%    delete(gca);

    dlm=strfind(sfiles(i).name,'.');
    fnameshort = sfiles(i).name(1:(dlm-1));
    fnamefull = sfiles(i).name(1:(dlm(end)-1));

    ffname = strcat(indir,filesep,sfiles(i).name); 

    strcat('Loaded file: ',ffname)


    ms0=loadmolxyz(ffname,fnameshort);
    

    ms0 = createbondtable(ms0);
    [ms0,status] = identmol(ms0,moltype);
    if status
      lastwarn
      continue
    end
    ms0 = createzmt(ms0);
    ms0 = zmt2xyz(ms0);

%    plotmol(ms0);
%    pause


    if exist(odir)~=7
       mkdir(odir);
    end

%    [xxx,order]=sortrows(ms0.pind);
%    order=[1 2 3 4 6];
    savemol(odir,ms0,0);
    savemol(odir,ms0,2);

 
%    'Press any key.'
%    pause
end

