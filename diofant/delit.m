format compact
max_factors = 1;
n_max_factors = 2;
for n=2:999999
 n_factor=0;
 for i=2:int32(n/2)
   if ~mod(n,i)
     n_factor=n_factor+1;
   end     
 end
 if n_factor>max_factors 
    max_factors = n_factor
    n_max_factors = n
 end
end
disp(n_max_factors);