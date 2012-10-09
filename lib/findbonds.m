function indvect = findbonds(mol,ind,marked)
%returns indexes of all bond of atom with index ind in molecule mol using
%bond table mol.bt
% marked is vector with indexes of atoms already proceeded
%
% Version 1.0    
% Last modified  R O Zhurakivsky 2007-08-04 
% Created        R O Zhurakivsky 2005-09-?

if nargin<3
    marked=[];
end


indvect=[];
for i=ind
  indvect=[indvect mol.btB(mol.btA==i) mol.btA(mol.btB==i)];
end

if ~isempty(marked)
    indvect=intersect(marked,uint16(indvect));
end

[xxx,I]=sort(mol.pind(indvect)); %#ok
indvect=indvect(I);

