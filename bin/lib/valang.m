function ang = valang(mols,ind1,ind2,ind3)
%returns planar angle between three points in the plane
%
% Version 1.0    
% Last modified  R O Zhurakivsky 2005-11-03
% Created        R O Zhurakivsky 2005-09-?

%when points are specified in different molecule structures then mols is 
%structure array
%otherwise it is sipmly molecule structure and ind1,ind2,ind3 are hardindexes 
%in it

%20050805 changed v1 vector orientation and acos argument sign


if numel(mols)~=1
    mol1=mols(1); mol2=mols(2); mol3=mols(3);
else
    mol1=mols; mol2=mols; mol3=mols;
end;

i1=find(mol1.ind==ind1); %true indexes
i2=find(mol2.ind==ind2);
i3=find(mol3.ind==ind3);
  
v1 = [-mol2.x(i2)+mol1.x(i1) -mol2.y(i2)+mol1.y(i1) -mol2.z(i2)+mol1.z(i1)];
v2 = [mol3.x(i3)-mol2.x(i2) mol3.y(i3)-mol2.y(i2) mol3.z(i3)-mol2.z(i2)];
ang = acos(dot(v1,v2)/norm(v1)/norm(v2))*180/pi;
