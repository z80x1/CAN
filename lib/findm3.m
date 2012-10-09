format compact;
i0=33;
for j=1:inf
  i = m3(i0);
  disp(i);
  if i==i0, break, end;
  i0 = i;
end    