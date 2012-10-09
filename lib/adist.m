function d = adist(mols,ind1,ind2)
%returns distance between two atoms
%
% Version 1.0    
% Last modified  R O Zhurakivsky 2005-09-20
% Created        R O Zhurakivsky 2005-09-20

%when points are specified in different molecule structures then mols is 
%structure array
%otherwise it is sipmly molecule structure and ind1,ind2 are hardindexes 
%in it

  if numel(mols)~=1
    d = norm(createvect(mols(1),ind1,mols(2),ind2));
  else
    d = norm(createvect(mols,ind1,mols,ind2));
  end;
