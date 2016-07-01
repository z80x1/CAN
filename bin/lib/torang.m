function ang = torang(mols,ind1,ind2,ind3,ind4)
%returns torsion angle created by four points (degrees)

%when points are specified in different molecule structures then mols is 
%structure array
%otherwise it is sipmly molecule structure and ind1,ind2,ind3,ind4 are hardindexes 
%in it
%
% Version 1.0    
% Last modified  R O Zhurakivsky 2005-11-13
% Created        R O Zhurakivsky 2005-09-?

if numel(mols)~=1
  mol1=mols(1); mol2=mols(2); mol3=mols(3); mol4=mols(4);
else
  mol1=mols; mol2=mols; mol3=mols; mol4=mols;
end;

i1=find(mol1.ind==ind1); %true indexes
i2=find(mol2.ind==ind2);
i3=find(mol3.ind==ind3);
i4=find(mol4.ind==ind4);

v1 = [mol2.x(i2)-mol1.x(i1) mol2.y(i2)-mol1.y(i1) mol2.z(i2)-mol1.z(i1)];
v2 = [mol3.x(i3)-mol2.x(i2) mol3.y(i3)-mol2.y(i2) mol3.z(i3)-mol2.z(i2)];
v3 = [mol4.x(i4)-mol3.x(i3) mol4.y(i4)-mol3.y(i3) mol4.z(i4)-mol3.z(i3)];

u1 = cross(v1,v2);
u2 = cross(v2,v3);

ang = acos(dot(u1,u2)/norm(u1)/norm(u2))*180/pi;


if det([u1; u2; v2])<0
  ang = -ang;
end
