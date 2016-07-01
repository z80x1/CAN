function mol=zmt2xyz(mol,order)
%converts internal coordinates (Z-matrix) to Cartesian ones
%needed mol fields:
%iR, R, ialfa, alfa, ibeta, beta
%order vector with hard indexes for Zmatrix creation order
%
% Version 1.0    
% Last modified  R O Zhurakivsky 2006-01-28
% Created        R O Zhurakivsky 2006-?-?
    
if nargin<2
   order=[];
%   [xxx,order]=sortrows(mol.pind); 
end
%order=reshape(order,1,numel(order));
or=createbondchain(mol,order);

pifac = pi/180;

mol.x(or(1))=0; mol.y(or(1))=0; mol.z(or(1))=0;

mol.x(or(2))=0; mol.y(or(2))=0; mol.z(or(2))=mol.R(or(2));

mol.x(or(3))=mol.R(or(3))*sin(mol.alfa(or(3))*pifac);
mol.y(or(3))=0;
mol.z(or(3))=mol.R(or(2))-mol.R(or(3))*cos(mol.alfa(or(3))*pifac);

for i=or(4:end)
	%from NWChem4.6 geom/geom_hnd.F

    iC = or(mol.iR(i));
    iB = or(mol.ialfa(i));
    iA = or(mol.ibeta(i));
    CB = [mol.x(iB)-mol.x(iC) mol.y(iB)-mol.y(iC) mol.z(iB)-mol.z(iC)];
    CB=CB/norm(CB);
    CA = [mol.x(iA)-mol.x(iC) mol.y(iA)-mol.y(iC) mol.z(iA)-mol.z(iC)];
    CA=CA/norm(CA);
    DOT = dot(CB,CA);

    CA=CA-DOT*CB;
    CA=CA/norm(CA);

    T=[CB; CA; cross(CB,CA)]';

    RCD=mol.R(i); aalfa=mol.alfa(i)*pifac; bbeta=-mol.beta(i)*pifac;
    xxD=RCD*cos(aalfa);
    yyD=RCD*sin(aalfa)*cos(bbeta);
    zzD=RCD*sin(aalfa)*sin(bbeta);

    mol.x(i)=mol.x(iC)+sum(T(1,1)*xxD+T(1,2)*yyD+T(1,3)*zzD);
    mol.y(i)=mol.y(iC)+sum(T(2,1)*xxD+T(2,2)*yyD+T(2,3)*zzD);
    mol.z(i)=mol.z(iC)+sum(T(3,1)*xxD+T(3,2)*yyD+T(3,3)*zzD);


end

