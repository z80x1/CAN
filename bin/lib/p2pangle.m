function ang = p2pangle(N1,N2)
%returns angle between planes with vectors N1 and N2 (in radians)
%
% Version 1.0    
% Last modified  R O Zhurakivsky 2006-02-25 
% Created        R O Zhurakivsky 2005-09-?

ang=-acos(dot(N1,N2)/norm(N1)/norm(N2));

%sgn=sign(prod(cross(N1,N2)))
