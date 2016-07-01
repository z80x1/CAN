function [i1,i2,i3]=find4alink(mol,i)
%returns indexes of 4 atoms that forms link starting from i and are marked
%051105: removed transpose mark from i1=ii1 and i2=ii2
%090419: if possible don't use links with glycosydic bond
%
% Version 1.0    
% Last modified  R O Zhurakivsky 2007-10-07 
% Created        R O Zhurakivsky 2005-09-?

pindsdef;
atomsind;

marked = find(mol.marked)';

variants=[];

ii1 = findbonds(mol,i,marked);
for i1=ii1

  marked=setdiff(marked,i1);
  ii2= findbonds(mol,i1,marked);
  for i2=ii2

    marked=setdiff(marked,i2);       
    ii3 = findbonds(mol,i2,marked);

    if numel(ii3)>0
      for i3=ii3
        
          variants(end+1,:)=[i i1 i2 uint16(i3)];
      end
    end 

  end
end

V=reshape(mol.pind(variants),size(variants));
A=sort(V,2);
B=sort(magictorsion,2);
[C,I]=intersect(A,B,'rows');
if isempty(I)
    
    pindpC1=find(strcmp(pind.labels,'pC1'));
    pindbN1=find(strcmp(pind.labels,'bN1'));
    pindbN9=find(strcmp(pind.labels,'bN9'));
    Vtor=sum((V==pindpC1)|(V==pindbN1)|(V==pindbN9),2);% find links with glicosydic bond atoms C1',N1 or N9
    if (sum(Vtor==2)<size(Vtor,1)) %if there are links without glicosydic bond - use them, just remove rest of links
        variants=variants(find(Vtor~=2),:);
    end    
    
    I=1; 
end

i1=variants(I,2);
i2=variants(I,3);
i3=variants(I,4);


