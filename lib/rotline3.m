function mol2 = rotline3(vect1,mol1,fi,aindex)
%ROTLINE3 rotates system of point around vector vect1 by angle fi 
%   vect1 - start vector of rotation
%   mol1.x,mol1.y,mol1.z - coordinates of system of point to rotate
%   fi - angle to rotate
%   aindex - index atom in vactors cx1,mol1.y,mol1.z which is the center of rotation
%
% Version 1.0    
% Last modified  R O Zhurakivsky 2006-08-20
% Created        R O Zhurakivsky 2006-?-?

aind=find(mol1.ind==aindex); %true index
mol2 = mol1;

j3 = vect1/norm(vect1);
j2 = cross(j3,[1 1 1]);
j2 = j2/norm(j2);
j1 = cross(j2,j3);

%movement of cytidine coordinate system to make aind atom be (0,0,0)
cx2 = mol1.x-mol1.x(aind);
cy2 = mol1.y-mol1.y(aind);
cz2 = mol1.z-mol1.z(aind);

cx3 = cx2*j1(1)+cy2*j1(2)+cz2*j1(3);
cy3 = cx2*j2(1)+cy2*j2(2)+cz2*j2(3);
cz3 = cx2*j3(1)+cy2*j3(2)+cz2*j3(3);
                      
cx4=cx3*cos(fi)-cy3*sin(fi);
cy4=cx3*sin(fi)+cy3*cos(fi);
cz4=cz3;

cx5 = cx4*j1(1)+cy4*j2(1)+cz4*j3(1);
cy5 = cx4*j1(2)+cy4*j2(2)+cz4*j3(2);
cz5 = cx4*j1(3)+cy4*j2(3)+cz4*j3(3);

mol2.x = cx5+mol1.x(aind);
mol2.y = cy5+mol1.y(aind);
mol2.z = cz5+mol1.z(aind);
