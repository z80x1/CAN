function val = strcmpcellar(cellar,str,flexact)
%finds strings equal to str in the string array cellar
%2006-0613: empty str detection added 
%
% Version 1.0    
% Last modified  R O Zhurakivsky 2007-08-24
% Created        R O Zhurakivsky 2006-?-?

if nargin<3
  flexact=1;
end


val=[];

if ~isempty(str)
  if flexact
      for i=1:numel(cellar)
        if strcmp(cellar{i},str)
          val(end+1)=i;
        end
      end
  else
      for i=1:numel(cellar)
        if strfind(cellar{i},str)
          val(end+1)=i;
        end
      end
  end
else
  for i=1:numel(cellar)
    if isempty(cellar{i})
      val(end+1)=i;
    end
  end
end
