%
% Version 1.0    
% Last modified  R O Zhurakivsky 2007-11-09
% Created        R O Zhurakivsky 2006-?-?

indir='D:\_diplom\CAN\xyz\rx19\rd';
moltype=13;


sfiles = dir(strcat(indir,filesep,'*.xyz'));
numfiles = size(sfiles,1);
if ~numfiles, error('No files found!'), end

for i=1:numfiles

    fname=sfiles(i).name;
    dlm=strfind(fname,'.');
    fnameshort = fname(1:(dlm-1));


    maxind=0;
    ms0=loadmolxyz([indir filesep fname],fnameshort,maxind);
    ms0 = createbondtable(ms0);
    [ms0,status] = identmol(ms0,moltype);
    ms0 = createzmt(ms0);

    orderanchor=[];
    %orderanchor=[4 6 1]; %pinds pC4, pO4, pC1
    order=[];
    for I=1:numel(orderanchor)
        order(end+1)=find(ms0.pind==orderanchor(I));
    end

    ms0.desc=[fnameshort '_out'];
    savemol(indir,ms0,0,order);
    savemol(indir,ms0,3,order);

end