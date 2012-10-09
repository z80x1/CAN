function [iprev,inext,sdesc_prev,sdesc_next] = checksimdesc(descstr,desc)

%CHECK IF CONFORMATION identifier IS NEW ONE
%desc is list of identifiers of existing conformations
%returns 
%iprev: number of conformation in list desc the descstr is similar to in previous P angle family or 0
%inext: number of conformation in list desc the descstr is similar to in next P angle family or 0
%sdesc_prev: string identifier of conformation with index iprev
%sdesc_next: string identifier of conformation with index inext
%
% Version 1.0    
% Last modified  R O Zhurakivsky 2007-04-07 
% Created        R O Zhurakivsky 2005-09-?

global GLsdescr
atomsind

sdesc_prev=descstr;
sdesc_prev(1)=sdesc_prev(1)-1;
if sdesc_prev(1)<GLsdescr(1)
    sdesc_prev(1)=GLsdescr(end);
end
%iprev=strcmpar(desc,sdesc_prev);
[c,iprev,ib]=intersect(desc,sdesc_prev,'rows');
if isempty(iprev), iprev=0;, end;

sdesc_next=descstr; 
sdesc_next(1)=sdesc_next(1)+1;
if sdesc_next(1)>GLsdescr(end)
    sdesc_next(1)=GLsdescr(1);
end
%inext=strcmpar(desc,sdesc_next);
[c,inext,ib]=intersect(desc,sdesc_next,'rows');
if isempty(inext), inext=0;, end;


