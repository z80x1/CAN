function y = srd(x,varargin)
%MEAN   Average or mean value.
%   For cells Y = STD(X) returns the standard deviation of the elements in each of cells
%
% Version 1.0    
% Last modified  R O Zhurakivsky 2007-11-29
% Created        R O Zhurakivsky 2007-11-28

y=nan(numel(x),1);
for i=1:numel(x)
    res=std(x{i},varargin{:});
    if ~isempty(res)
        y(i)=res;
    else
        y(i)=nan;
    end
end
y=reshape(y,size(x));

