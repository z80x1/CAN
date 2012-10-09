function [ind2,dist]=findnearestatom(mol,atomtype,ind1)
%find the nearest atom of type atomtype to ind1 atom
%return index of finded atom and distance
%
% Version 1.0    
% Last modified  R O Zhurakivsky 2006-07-08 
% Created        R O Zhurakivsky 2005-09-?

iatoms = strcmpcellar(mol.labels,atomtype);
if isempty(iatoms)
  error('findnearestatom: no required atoms found2')
end

distvect=[mol.x(iatoms)-mol.x(ind1) mol.y(iatoms)-mol.y(ind1) mol.z(iatoms)-mol.z(ind1)];


dists=sqrt(distvect(:,1).^2+distvect(:,2).^2+distvect(:,3).^2);

[dist,indX]=min(dists);

ind2=iatoms(indX);

if ~dist
  error('findnearestatom: zero distance detected')
end