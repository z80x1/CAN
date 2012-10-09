function [ms1,errlev] = rotaroundbond(ms0,ind,angle)
%
% Version 1.0    
% Last modified  R O Zhurakivsky 2007-08-20
% Created        R O Zhurakivsky 2006-?-?

%rotates atoms with indexes in ind by value 'angle'

ms1=ms0;
errlev=1;

N=numel(ind);
if ~N
  warning('rotaround:emptyvector','Empty vector with indexes of atoms to rotate');
  return;
end
ind = reshape(ind,[1,N]);

for i=ind(1:N)
   ms1.beta(i) = mod(ms1.beta(i)+angle,360);
end
ms1=zmt2xyz(ms1);
errlev=0;