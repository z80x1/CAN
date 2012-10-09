function    drawbond(mol1,ind1,mol2,ind2)
%draws bond between atoms with indexes ind1 and ind2
%
% Version 1.0    
% Last modified  R O Zhurakivsky 2005-09-20 
% Created        R O Zhurakivsky 2005-09-?

global flPlot

if flPlot~=0
    return
end


i1=find(mol1.ind==ind1); %true indexes
i2=find(mol2.ind==ind2);

a1 = [mol1.x(i1) mol1.y(i1) mol1.z(i1)];
a2 = [mol2.x(i2) mol2.y(i2) mol2.z(i2)];

plot3([a1(1) a2(1)],[a1(2) a2(2)],[a1(3) a2(3)],'k-');
grid on
rotate3d on
