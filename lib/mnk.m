function [a, b, err] = mnk(x,f)
%MNK solve linear Minimal Quatrature Method
% f(x) = a + bx
% err - rmse 


n = numel(x);

xx = mean(x);
ff = mean(f);


b = sum((x-xx).*(f-ff))/sum((x-xx).^2);

a = ff - b*xx;

err=rmse(f, b*x+a);