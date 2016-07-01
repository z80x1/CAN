function [x,y,z]=sfer2dec(r,teta,fi)
% transform spherical coordinates to cartesian ones
%
% Version 1.0    
% Last modified  R O Zhurakivsky 2005-11-03
% Created        R O Zhurakivsky 2005-09-?

%teta=mod(teta,pi);
x=r.*sin(teta).*cos(fi);
y=r.*sin(teta).*sin(fi);
z=r.*cos(teta);
