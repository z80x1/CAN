function dist=dist2linecart(P0,P1,P2)
%returns distance from point P0 to line P1P2
%
% Version 1.0    
% Last modified  R O Zhurakivsky 2007-08-19 
% Created        R O Zhurakivsky 2005-09-?

%   d(P,L) = |uL X w|
%   uL = vL / |vL|
%   vL = P1P2,  w = P1P0

vL = P2-P1;
uL = vL/norm(vL);

%disp([uL P0-P1])
dist = norm(cross(uL, P0-P1));

