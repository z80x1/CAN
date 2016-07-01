function ms0=loadmolxyz(ffname,desc,maxind)
%load molecule geometry from xyz file
%
% Version 1.0    
% Last modified  R O Zhurakivsky 2007-08-04 
% Created        R O Zhurakivsky 2005-09-?

if nargin<3
    maxind=0;
end
if nargin<2
    desc='NULL desc';
end

[ms0.labels,ms0.x,ms0.y,ms0.z]=textread(ffname,'%s%f%f%f');
ms0.labels=cell(ms0.labels);
%ms0.x=x';
%ms0.y=y';
%ms0.z=z';
ms0.atomnum = uint16(length(ms0.labels));
ms0.ind = ((maxind+1):(maxind+ms0.atomnum))';
ms0.desc = desc;
