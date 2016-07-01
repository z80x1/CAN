function ind = identatom(mol,atype,ilinks,tlinks,tnolinks)
%identatom(ms0,'O',[iC5])
%returns index of specified atom of type atype in molecule mol using bond table 
%bt that have links with atoms with indexes specified in ilinks
% and have links with atoms of types specified in tlinks
% and haven't links with atoms of types specified in tnolinks

% by now only one element in ilinks and tlinks is supported
%if failed returns 0
%
% Version 1.0    
% Last modified  R O Zhurakivsky 2007-02-25
% Created        R O Zhurakivsky 2005-09-?

if nargin<5
    tnolinks='';
end
if nargin<4
    tlinks='';
end


ind = findatom(mol,atype,1);
while ind~=0
  idxs = findbonds(mol,ind);
  if sum(idxs==ilinks(1))==1
    if isempty(tlinks) || numel(strcmpcellar(mol.labels(idxs),tlinks)) > 0 % check if found atom have links with atoms of type tlinks
      if isempty(tnolinks) || numel(strcmpcellar(mol.labels(idxs),tnolinks)) == 0 % check if found atom have links with atoms of type tlinks
        break
      end
    end
  end
  ind = findatom(mol,atype,ind+1);
end
