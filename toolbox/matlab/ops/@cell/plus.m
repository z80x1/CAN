function C = PLUS(A,B)
%+   Plus.
%   X + Y adds matrices X and Y.  X and Y must have the same
%   dimensions unless one is a scalar (a 1-by-1 matrix).
%   A scalar can be added to anything.  
%
%   C = PLUS(A,B) is called for the syntax 'A + B' when A or B is an
%   object.
%MIN    Smallest component.
%   For cells returns the smallest elements for each cell
%
% Version 1.0    
% Last modified  R O Zhurakivsky 2007-11-29
% Created        R O Zhurakivsky 2007-11-29

if iscell(B)
  error('Dont know what to do with two cell objects');
  return
end

for i=1:numel(A)
    C{i}=A{i}+B;
end
