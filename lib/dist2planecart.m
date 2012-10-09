function dist=dist2planecart(P0,P1,P2,P3)
%returns distance from point with hardindex a0 in molecule structure mol
%to plane than consist points with indexes a1, a2, a3)
%
% Version 1.0    
% Last modified  R O Zhurakivsky 2007-08-10 
% Created        R O Zhurakivsky 2005-09-?

%    coefA = at1y*(at2z-at3z)+at2y*(at3z-at1z)+at3y*(at1z-at2z)
%    coefB = at1z*(at2x-at3x)+at2z*(at3x-at1x)+at3z*(at1x-at2x)
%    coefC = at1x*(at2y-at3y)+at2x*(at3y-at1y)+at3x*(at1y-at2y)
%    prop.dplC2=(coefA*(at0x-at1x)+coefB*(at0y-at1y)+coefC*(at0z-at1z))/sqrt(coefA*coefA+coefB*coefB+coefC*coefC)

%2007-0810 plot distance from point P1 to plane which consist points P2,P3,P4. Px are 3D vectors

T1=[P1(2) P2(2) P3(2); ...
    P1(3) P2(3) P3(3); ...
    P1(1) P2(1) P3(1)];

T2=[P2(3) P3(3) P1(3); ...
    P2(1) P3(1) P1(1); ...
    P2(2) P3(2) P1(2)];

T3=[P3(3) P1(3) P2(3); ...
    P3(1) P1(1) P2(1); ...
    P3(2) P1(2) P2(2)];

plcoef = sum(T1.*(T2-T3),2);

a1 = [P1(1) P1(2) P1(3)];
a0 = [P0(1) P0(2) P0(3)];
dist = dot(plcoef,a0-a1)/norm(plcoef);
