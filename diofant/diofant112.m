format compact
n=1;
for i=1:2010
  n=rem(n*9,1000);
%  disp([int2str(i)  ' : '  int2str(n)]);
end
n