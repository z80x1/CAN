function vect = createvect(mol1,ind1,mol2,ind2)
%  returns the distance vector between atoms with indexes ind1 and ind2 in 
%molecules mol1 & mol2 accordingly
%
% Version 1.0    
% Last modified  R O Zhurakivsky 2006-08-20 
% Created        R O Zhurakivsky 2005-09-?

i1=find(mol1.ind==ind1); %true indexes
i2=find(mol2.ind==ind2);
if isempty(i1) || isempty(i2)
  error('indexes not found');
end

vect = [mol2.x(i2)-mol1.x(i1) mol2.y(i2)-mol1.y(i1) mol2.z(i2)-mol1.z(i1)];
