function val = strcmpar(ar,str)
%finds string equal to str in array ar
%
% Version 1.0    
% Last modified  R O Zhurakivsky 2005-09-20
% Created        R O Zhurakivsky 2005-09-?

val=0;
for i=1:size(ar,1)
  if strcmp(ar(i,:),str)
    val=i;
    return
  end
end