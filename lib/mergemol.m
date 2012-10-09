function m3=mergemol(m1,m2)
%merge two molecules
%
% Version 1.0    
% Last modified  R O Zhurakivsky 2006-10-04 
% Created        R O Zhurakivsky 2006-?-?

m3.atomnum = uint16(m1.atomnum) + uint16(m2.atomnum);
m3.labels = [m1.labels ; m2.labels];
%m3.ind = [m1.ind ; m2.ind];
m3.ind = 1:m3.atomnum;


m3.x = [m1.x ; m2.x];
m3.y = [m1.y ; m2.y];
m3.z = [m1.z ; m2.z];
