function [r,teta,fi]=dec2sfer(x,y,z)
% transform Decart coordinates to spherical ones
%
% Version 1.0    
% Last modified  R O Zhurakivsky 2006-08-20 
% Created        R O Zhurakivsky 2005-09-?

r=sqrt(x.^2+y.^2+z.^2);

warning off MATLAB:divideByZero
teta=acos(z./r);
teta(isnan(teta))=0;

fi=atan(y./x)+(x<0)*pi;
fi(isnan(fi))=0;
warning on MATLAB:divideByZero
