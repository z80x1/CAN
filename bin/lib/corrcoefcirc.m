function rc = corrcoefcirc(X)
%calculate circular correlation coefficient matrix accordingly to 
%Kitamura, K., Wakahara, A., Mizuno, H., Baba, Y. & Tomita, K. (1981). Conformationally ``concerted'' 
%changes in nucleotide structures. A new description using circular correlation and regression analyses. 
%J. Am. Chem. Soc. 103, 3899-3904.
%assumes angles are in radians
%
% Version 1.0    
% Last modified  R O Zhurakivsky 2006-06-13 
% Created        R O Zhurakivsky 2005-09-?


n=size(X);

alfa = sum(cos(X))/n(2);
beta = sum(sin(X))/n(2);
omega = sqrt(alfa.^2+beta.^2);
theta0 = atan(beta./alfa);
theta00 = repmat(theta0,n(1),1);

X = X - theta00;

rc=zeros(n(2));
for i=1:n(2)
for j=i:n(2)
    rc(i,j) = sum(sin(X(:,i)).*sin(X(:,j)))/sqrt(sum(sin(X(:,i)).^2)*sum(sin(X(:,j)).^2));
end
end

rc = rc + rc' - diag(diag(rc));