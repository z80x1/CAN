function  atomchain = createbondchain(mol,order,fl_Hlast)
%builds chain of atoms where any following atom can be characterized by previous ones.
%returns physical order in mol structure
%fl_Hlast=1 - will put indexes of H atoms at the end
%
% Version 1.0    
% Last modified  R O Zhurakivsky 2007-08-07 
% Created        R O Zhurakivsky 2005-09-?

if nargin<1
    error('no molecular structure specified');
end
if nargin < 2 || isempty(order) || ~all(order)
  order=1; %if no initial order specified start with the first atom
end
if nargin < 3
  fl_Hlast=0;
end

atomchain = uint16(zeros(1,mol.atomnum));

fillednum=numel(order);
atomchain(1:fillednum) = order;

flbreak=0; %flag to break after numeration disorder occurs
j=0;

fl_allheavydetected=0;
Hinds=find(ismember(mol.labels,'H'));
Nainds=find(ismember(mol.labels,'Na'));
num_heavyatoms=mol.atomnum-numel(Hinds)-numel(Nainds);

%if ~isfield(mol,'btA') || ~isfield(mol,'btB')
%    mol.btA=zeros(0,0,'uint16');
%    mol.btB=zeros(0,0,'uint16');
%end

while 1
  restind=setdiff(1:mol.atomnum,atomchain);
  if isempty(restind)
     break
  end

  if ~fl_allheavydetected && numel(find(atomchain))>=num_heavyatoms %zhr100218
      fl_allheavydetected=1;
  end
  if fl_Hlast && ~fl_allheavydetected
      %skip if not all heavy atoms included in chain yet
      restind=setdiff(restind,Hinds);
  end

 % if isfield(mol,'pind')
      [xxx,I]=sort(mol.pind(restind));
      inds = restind(I);
  %else
  %    inds=restind;
  %end

  for i=inds

%    ii = intersect(atomchain,findbonds(mol,i));
    ii = findbonds(mol,i,atomchain);
    if ~isempty(ii)
      
      atomchain(fillednum+1)=i;
      fillednum = fillednum+1; 
      flbreak=0;
      break
    else
      flbreak=1;  
    end
  end

  j=j+1;
  if j>1e3
    error('createbondchain: infinite loop detected');
  end    
end
return