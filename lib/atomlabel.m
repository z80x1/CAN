function str=atomlabel(mol,i)
%returns label for atom i in mol structure
%
% Version 1.0    
% Last modified  R O Zhurakivsky 2012-04-23 
% Created        R O Zhurakivsky 2012-04-23

str = [mol.labels{i} int2str(i)];

