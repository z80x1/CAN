function mol=createbondtable(mol,mode)
%creates table of chemical bonds from molecule geometry mol
%
%mode
%if 1 - detect Hbonds - distancies between hydrogen and other atom less
%than 1.3A
% Version 1.0    
% Last modified  R O Zhurakivsky 2007-08-04 
% Created        R O Zhurakivsky 2005-09-?

if nargin<2
    mode=0;
end


A=zeros(0,0,'uint16');
B=zeros(0,0,'uint16');

for i=1:mol.atomnum
  if strcmp(mol.labels(i),'Na')
      fl_iNa=1;
  else
      fl_iNa=0;
  end;
  if strcmp(mol.labels(i),'P')
      fl_iP=1;
  else
      fl_iP=0;
  end
  
  for j=(i+1):mol.atomnum
      
     %create fictive bond (only one) for Na atom
     if (fl_iNa || strcmp(mol.labels(j),'Na')) 
         fl_Na = 1;
     else
         fl_Na = 0;
     end     
     if (fl_iP || strcmp(mol.labels(j),'P')) 
         fl_P = 1;
     else
         fl_P = 0;
     end     
     
    if (fl_Na && fl_iP) || ...
            (~fl_Na && adist(mol,i,j)<=bondlen(mol.labels(i),mol.labels(j))) || ...
            (mode==1 && adist(mol,i,j)<=1.5) %
      A(end+1) = i;
      B(end+1) = j;
    end
  end
end

mol.btA=A;
mol.btB=B;