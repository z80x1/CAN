%
% Version 1.0    
% Last modified  R O Zhurakivsky 2006-03-10
% Created        R O Zhurakivsky 2006-03-10

atomsind

x=[
3817	%O3'H       4
3684    %O5'H...O2  2
3836    %O5'H       5
3766	%O3'H       3
3842    %O5'H       6
3617 %N3H           1
];

y=[
0.3
0.6
0.5
0.1
0.55
1
];


x=x*CC.freqcoefB3LYP631Gdp;
xx=[x x]';

yy=[zeros(size(y)) y]';

plot(xx,yy,'k')
dx=max(x)-min(x);






N=5000;
sumY=zeros(1,N);

  xlen=max(x)-min(x);
minx=min(x)-0.1*xlen;
maxx=max(x)+0.1*xlen;





  X=minx:(maxx-minx)/(N-1):maxx;



  for I=1:N
    Y(I)=sum(y.*exp(-0.07*(X(I)-x).^2));
%    Ynorm(I)=sum(ynorm.*exp(-0.1*(X(I)-x).^2));

  end
%  sumY=sumY+Ynorm;


plot(X,Y)

axis([min(X) max(X) 0 1.1*max(Y)])