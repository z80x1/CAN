function    vect = createplane(mol1,ind1,mol2,ind2,mol3,ind3)
%returns vector of the plane where atoms with indexes ind1, ind2 and ind3 from
%molecules mol1, mol2  and mol3 lays
%
% Version 1.0    
% Last modified  R O Zhurakivsky 2006-08-20 
% Created        R O Zhurakivsky 2005-09-?

i1=find(mol1.ind==ind1); %true indexes
i2=find(mol2.ind==ind2);
i3=find(mol3.ind==ind3);

A = det([mol2.y(i2)-mol1.y(i1) mol2.z(i2)-mol1.z(i1);...
		 mol3.y(i3)-mol1.y(i1) mol3.z(i3)-mol1.z(i1)]);
B = det([mol2.z(i2)-mol1.z(i1) mol2.x(i2)-mol1.x(i1);...
		 mol3.z(i3)-mol1.z(i1) mol3.x(i3)-mol1.x(i1)]);
C = det([mol2.x(i2)-mol1.x(i1) mol2.y(i2)-mol1.y(i1);...
		 mol3.x(i3)-mol1.x(i1) mol3.y(i3)-mol1.y(i1)]);
vect = [A B C];
