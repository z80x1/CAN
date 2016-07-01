function h1g = plotspectra(x,y,color)
%plotspectra plot spectra from modes
% x - frequencies vector
% y - intencities vector
%
% Version 1.0    
% Last modified  R O Zhurakivsky 2009-04-12
% Created        R O Zhurakivsky 2006-04-12

if nargin < 2
    error('Need 2 vector with data');
end
if nargin < 3
    color='r';
end


N=4000; %number of points to calculate spectra in
halfwidth=8.5; %cm^-1 - halfwidth of gaussian
a=4*log(2)/halfwidth^2;
xlen=1.15*max(x);
X=0:xlen/(N-1):xlen;
Y=zeros(size(X));
for I=1:N
    Y(1,I)    =sum(y.*    exp(-a*(X(I)-x).^2));
end
h1g=plot(X,Y(1,:),color);
%set(h1g,'LineWidth',2);
