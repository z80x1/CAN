%plot transition diagram 
%
% Version 1.0    
% Last modified  R O Zhurakivsky 2006-07-13
% Created        R O Zhurakivsky 2006-?-?

atomsind

y=[2.31 2.81 1.55 2.58 3.48 2.92] %#ok %MP2_6311__Gdp (TEC corrected) 
eps=[-65.8 60.1 174.6 -120.9 -1.9 111.5] %#ok

eps0=[-65.8 -120.9 60.1 -1.9 111.5 174.6] %#ok
y0=[2.84 3.99 2.95 4.48 3.87 2.27] %#ok % DFT_631Gdp

[eps,ind]=sort(eps);
[eps0,ind0]=sort(eps0);
y=y(ind);
eps(end+1)=eps(1)+360;
y(end+1)=y(1);

y0=y0(ind0);
eps0(end+1)=eps0(1)+360;
y0(end+1)=y0(1);


plot([eps;eps0]',[y;y0]')
axis([min(eps) max(eps) 0 1.05*max([y y0])])
xlabel('eps,degree')
ylabel('dE_{6-31Gdp},kcal/mol')

legend('MP2_{6311++Gdp (TEC corrected)}','DFT_{631Gdp}');

