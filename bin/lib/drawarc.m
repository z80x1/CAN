function drawarc(a,txt)
%plots arcs on conformational ring plot 
%a  - [a1 a2] - arc limits in radians
%txt - arc label
%
% Version 1.0    
% Last modified  R O Zhurakivsky 2007-08-04 
% Created        R O Zhurakivsky 2005-09-?

R=1.25; %arc radius
dR=0.03; %halflength if limiting segment
dRtxt=0.10;
dRshift = 0.03;

a1=min(a);
a2=max(a);
tfi=a1:pi/400:a2;
tr=repmat(R,size(tfi));
%tr=tr+((mod(tfi,2*pi)>0.5*pi)*(mod(tfi,2*pi)<1.5*pi)*2-1)*dRshift;

plot(tr.*cos(-tfi+pi/2),tr.*sin(-tfi+pi/2),'k'); %draws arc itself

fi=-tfi(1)+pi/2; 
X=[(R-dR).*cos(fi); (R+dR).*cos(fi)]; %first arc limiting segment
Y=[(R-dR).*sin(fi); (R+dR).*sin(fi)]; 
plot(X,Y,'k');

fi=-tfi(end)+pi/2;
X=[(R-dR).*cos(fi); (R+dR).*cos(fi)]; %second arc limiting segment
Y=[(R-dR).*sin(fi); (R+dR).*sin(fi)];
plot(X,Y,'k');

%070325
%fi=-tfi(int8(14*end/16))+pi/2;
fi=-tfi(int16(end/2))+pi/2;

Rtxt = R+dRtxt+((mod(fi,2*pi)>0.5*pi)*(mod(fi,2*pi)<1.5*pi)*2-1)*dRshift;
%Rtxt = R+dRtxt;
X=Rtxt.*cos(fi);
Y=Rtxt.*sin(fi);
text(X,Y,txt,'FontSize',36);

