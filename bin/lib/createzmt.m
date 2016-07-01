function mol = createzmt(mol,order)
%create Z-matrix fields in mol    
%order - hardindexes
%
% Version 1.0    
% Last modified  R O Zhurakivsky 2007-08-04 
% Created        R O Zhurakivsky 2005-09-?

%060128 added creation of bonds chain

if nargin<2
%   [xxx,order]=sortrows(mol.pind); 
  order=[];
end
%order=reshape(order,1,numel(order));
    
%or=[order,setdiff(1:mol.atomnum,order)];
or=createbondchain(mol,order);
%mol.pind(or)'

mol.zmtatomorder=or;

mol.marked=zeros(mol.atomnum,1,'uint16');
mol.iR=zeros(mol.atomnum,1,'uint16');
mol.ialfa=zeros(mol.atomnum,1,'uint16');
mol.ibeta=zeros(mol.atomnum,1,'uint16');
mol.R=zeros(mol.atomnum,1,'double');
mol.alfa=zeros(mol.atomnum,1,'double');
mol.beta=zeros(mol.atomnum,1,'double');

mol.iR(or(2))   = 1;
mol.R(or(2))    = adist(mol,or(2),or(1));

mol.iR(or(3))   = 2;
mol.ialfa(or(3)) = 1;
mol.R(or(3))    = adist(mol,or(3),or(2));
mol.alfa(or(3)) = valang(mol,or(3),or(2),or(1));

mol.marked(or(1:3))=1;
for i=or(4:end)
    [i1,i2,i3] = find4alink(mol,i);

    mol.iR(i) = find(or==i1);
    mol.R(i) = adist(mol,i,i1);
    mol.ialfa(i) = find(or==i2);
    mol.alfa(i) = valang(mol,i,i1,i2);
    mol.ibeta(i) = find(or==i3);
    mol.beta(i) = torang(mol,i,i1,i2,i3);
    mol.marked(i)=1;
end
