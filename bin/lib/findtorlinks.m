function I=findtorlinks(mol,i2,i3)
%returns indexes of 4 linked atoms that contain indexes i2 and i3
%for torsion rotation
%
% Version 1.0    
% Last modified  R O Zhurakivsky 2007-08-20 
% Created        R O Zhurakivsky 2005-09-?

ii1=setdiff(findbonds(mol,i2),i3);
ii4=setdiff(findbonds(mol,i3),i2);

[i1,i4]=meshgrid(ii1,ii4);
N=numel(i1);
i1=reshape(i1,N,1);
i4=reshape(i4,N,1);
i2=repmat(i2,N,1);
i3=repmat(i3,N,1);

I=[i1,uint16(i2),uint16(i3),i4];
