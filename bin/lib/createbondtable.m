function mol=createbondtable(mol)
%creates table of chemical bonds from molecule geometry mol
%
% Version 1.0    
% Last modified  R O Zhurakivsky 2007-08-04 
% Created        R O Zhurakivsky 2005-09-?

A=zeros(0,0,'uint16');
B=zeros(0,0,'uint16');

for i=1:mol.atomnum
  for j=(i+1):mol.atomnum
    if adist(mol,i,j)<=bondlen(mol.labels(i),mol.labels(j))
      A(end+1) = i;
      B(end+1) = j;
    end
  end
end

mol.btA=A;
mol.btB=B;