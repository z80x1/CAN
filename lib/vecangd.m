function ang = vecangd(mol,i1,i2,i3,i4)
%returns angle between vectors i1i2 and i3i4
%i1,i2,i3,i4 - hardindexes
%
% Version 1.0    
% Last modified  R O Zhurakivsky 2005-11-03
% Created        R O Zhurakivsky 2005-09-?

v1=[mol.x(i2)-mol.x(i1) mol.y(i2)-mol.y(i1) mol.z(i2)-mol.z(i1)];
v2=[mol.x(i4)-mol.x(i3) mol.y(i4)-mol.y(i3) mol.z(i4)-mol.z(i3)];
v3=[mol.x(i3)-mol.x(i1) mol.y(i3)-mol.y(i1) mol.z(i3)-mol.z(i1)];

ang=acosd(dot(v1,v2)/norm(v1)/norm(v2));

sign=dot(cross(v1,v2),v3);
if sign<0 
  ang=360-ang;
end
