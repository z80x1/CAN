function sdesc=rotatesdesc(sdesc,i)
%change character i description string sdesc on order a->b->c->a
%
% Version 1.0    
% Last modified  R O Zhurakivsky 2005-11-03
% Created        R O Zhurakivsky 2005-09-?

sdesc(i)=sdesc(i)+1;
if sdesc(i)>'c'
    sdesc(i)='a';
end
