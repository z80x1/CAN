function y = min(x,varargin)
%MIN    Smallest component.
%   For cells returns the smallest elements for each cell
%
% Version 1.0    
% Last modified  R O Zhurakivsky 2007-11-29
% Created        R O Zhurakivsky 2007-11-28

y=nan(numel(x),1);
for i=1:numel(x)
    res=min(x{i},varargin{:});
    if ~isempty(res)
        y(i)=res;
    else
        y(i)=nan;
    end
end
y=reshape(y,size(x));
