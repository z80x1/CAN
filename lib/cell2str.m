function str = cell2str(cellar,dlm)
%cell2str: compose string from cell array with dlm as delimiter
%
% Version 1.0    
% Last modified  R O Zhurakivsky 2011-10-27
% Created        R O Zhurakivsky 2011-10-27
if nargin<2
  dlm=' ';
end
if nargin<1
  error('no array specified');
end


str = '';
for i=1:numel(cellar)
    if i==1
        str = [str num2str(cellar{i})];
    else
        str = [str ' ' num2str(cellar{i})];
    end
    
end