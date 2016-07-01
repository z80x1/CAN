function ind = findatom(mol,atype,startind)
%returns index of first atom of type atype starting from startind
%
% Version 1.0    
% Last modified  R O Zhurakivsky 2009-04-19 
% Created        R O Zhurakivsky 2005-09-?

if nargin < 3
    error('Need 3 arguments');
end

for i=startind:mol.atomnum
%  if mol.labels{i}==atype && mol.marked(i)==0 %071029 change to comply
%  with multichar atom labels
  if strcmp(mol.labels(i),atype) && mol.marked(i)==0
    ind = i;
    return
  end
end
ind = 0;
