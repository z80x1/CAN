function V = remels(V,n)
%returns vector with removed elements with numbers provided in vector n
%
% Version 1.0    
% Last modified  R O Zhurakivsky 2005-09-20
% Created        R O Zhurakivsky 2005-09-20

N=numel(V);
V=reshape(V,N,1);

for i=1:numel(n)
  if n(i) == 1 
    V = V(2:N);
  elseif n(i) == N
    V = V(1:N-1);
  else
    V = [V(1:n(i)-1) ; V(n(i)+1:N)];
  end
  N=N-1;
  n=n-1;
end