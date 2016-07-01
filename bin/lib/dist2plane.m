function dist=dist2plane(mol,at0,at1,at2,at3)
%returns distance from point with hardindex a0 in molecule structure mol
%to plane than consist points with indexes a1, a2, a3)
%
% Version 1.0    
% Last modified  R O Zhurakivsky 2005-11-03 
% Created        R O Zhurakivsky 2005-09-?

%    coefA = at1y*(at2z-at3z)+at2y*(at3z-at1z)+at3y*(at1z-at2z)
%    coefB = at1z*(at2x-at3x)+at2z*(at3x-at1x)+at3z*(at1x-at2x)
%    coefC = at1x*(at2y-at3y)+at2x*(at3y-at1y)+at3x*(at1y-at2y)
%    prop.dplC2=(coefA*(at0x-at1x)+coefB*(at0y-at1y)+coefC*(at0z-at1z))/sqrt(coefA*coefA+coefB*coefB+coefC*coefC)

T1=[mol.y(at1) mol.y(at2) mol.y(at3); ...
    mol.z(at1) mol.z(at2) mol.z(at3); ...
    mol.x(at1) mol.x(at2) mol.x(at3)];

T2=[mol.z(at2) mol.z(at3) mol.z(at1); ...
    mol.x(at2) mol.x(at3) mol.x(at1); ...
    mol.y(at2) mol.y(at3) mol.y(at1)];

T3=[mol.z(at3) mol.z(at1) mol.z(at2); ...
    mol.x(at3) mol.x(at1) mol.x(at2); ...
    mol.y(at3) mol.y(at1) mol.y(at2)];
plcoef = sum(T1.*(T2-T3),2);

a1 = [mol.x(at1) mol.y(at1) mol.z(at1)];
a0 = [mol.x(at0) mol.y(at0) mol.z(at0)];
dist = dot(plcoef,a0-a1)/norm(plcoef);
