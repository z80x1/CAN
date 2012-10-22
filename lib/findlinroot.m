function [X0, ind] = findlinroot(X, Y)
%find one of the roots of discrete dependence Y(X)
%by linear interpolation
%X and Y are one dimension arrays
%X must be sorted 
%    
%returns:
%XO - interpolation result
%ind - index of point in X is nearest to 0
%
% Version 1.0    
% Last modified  R O Zhurakivsky 2012-10-15 
% Created        R O Zhurakivsky 2012-10-15

[m,i1] = min(abs(Y));
i1 = i1(1); %exclude duplicating points

if Y(i1)*Y(i1+1)<0
    i2 = i1 + 1;
else
    i2 = i1 - 1;
end

X0 = X(i1)-(X(i2)-X(i1))*Y(i1)/(Y(i2)-Y(i1));

ind = i1;

