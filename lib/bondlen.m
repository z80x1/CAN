function d = bondlen(a1,a2)
%returns typical bond length between atoms of types a1 & a2
%
% Version 1.0    
% Last modified  R O Zhurakivsky 2009-11-19 
% Created        R O Zhurakivsky 2005-09-?

global GLaspec

k=1.1;

d = k*(GLaspec.corad(strcmpcellar(GLaspec.type,a1))+GLaspec.corad(strcmpcellar(GLaspec.type,a2)));
