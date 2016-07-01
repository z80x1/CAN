function mol2 = delatom(mol1,index)
% deletes atom with index from molecule 
%
% Version 1.0    
% Last modified  R O Zhurakivsky 2007-08-04 
% Created        R O Zhurakivsky 2005-09-?

mol1.labels=reshape(mol1.labels,numel(mol1.labels),1);
ind=find(mol1.ind==index); %true index

mol2 = mol1;
if index == 1
    mol2.x = mol1.x(2:end);
    mol2.y = mol1.y(2:end);
    mol2.z = mol1.z(2:end);
    mol2.labels = mol1.labels(2:end);
    mol2.ind = mol1.ind(2:end);
    if isfield(mol1,'pind')
        mol2.pind = mol1.pind(2:end);
    end
    
elseif index == size(mol1.x)
    mol2.x = mol1.x(1:end-1);
    mol2.y = mol1.y(1:end-1);
    mol2.z = mol1.z(1:end-1);
    mol2.labels = mol1.labels(1:end-1);
    mol2.ind = mol1.ind(1:end-1);
    if isfield(mol1,'pind')
	  mol2.pind = mol1.pind(1:end-1);
    end
    
else
    mol2.x = [mol1.x(1:ind-1); mol1.x(ind+1:end)];
    mol2.y = [mol1.y(1:ind-1); mol1.y(ind+1:end)];
    mol2.z = [mol1.z(1:ind-1); mol1.z(ind+1:end)];
    mol2.labels = [mol1.labels(1:ind-1); mol1.labels(ind+1:end)];
    mol2.ind = [mol1.ind(1:ind-1); mol1.ind(ind+1:end)];
    if isfield(mol1,'pind')
	  mol2.pind = [mol1.pind(1:ind-1); mol1.pind(ind+1:end)];
    end
    
end

mol2.atomnum = uint16(length(mol2.labels));
