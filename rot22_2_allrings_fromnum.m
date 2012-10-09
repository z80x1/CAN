    %rot22_2: create "conformations rings" for beta, gamma, epsilon, chi and P angles at one picture
%
% Version 1.0    
% Last modified  R O Zhurakivsky 2008-07-31
% Created        R O Zhurakivsky 2006-?-?

atomsind
format compact

title='d4T'
ringtypes=['b'; 'g'; 'c'; 'p'];

tbeta=[-74.7
-73.4
80.0
89.6
-176.4
167.8
-62.0
-39.5
64.9
49.1
164.4
-53.7
-59.7
59.7
64.2
-170.5
-157.7
]';

tgamma=[
-65.3
-68.2
-68.7
-65.4
-66.1
-69.8
50.7
64.2
62.3
35.5
45.6
174.9
176.7
178.8
-178.8
-171.2
-171.0
]';

tchi=[
-109.5
67.9
-111.3
69.9
-109.3
68.1
-135.8
86.8
-107.6
64.9
-123.2
-110.2
70.6
-114.3
70.6
-113.3
68.9
]';

Pdeg=[
90.2
142.3
91.5
159.1
88.3
138.5
101.5
215.3
80.7
106.5
85.2
92.0
154.5
92.0
180.0
88.4
164.7
]';

cla
hold on

r1=0.2; %inner radius
dr=0.2; %outer radius
deltafi=30; %interval of fi labeling
tickdeltafi=30; %interval of tick placing
ticklen=0.02;

for ind=1:size(ringtypes)

    ringtype=ringtypes(ind);
    if strcmp(ringtype,'b')
      fi=tbeta*pi/180;
      mytitle='\beta';
    elseif strcmp(ringtype,'g')
      fi=tgamma*pi/180;
      mytitle='\gamma';
    elseif strcmp(ringtype,'d')
      fi=tdelta*pi/180;
      mytitle='\delta';
    elseif strcmp(ringtype,'e')
      fi=tepsilon*pi/180;
      mytitle='\epsilon';
    elseif strcmp(ringtype,'c')
      fi=tchi*pi/180;
      mytitle='\chi';
    elseif strcmp(ringtype,'p')
      fi=Pdeg*pi/180;
      %fi=P;
      mytitle='P';
    elseif strcmp(ringtype,'t')
      fi=tteta*pi/180;
      mytitle='\theta';
    elseif strcmp(ringtype,'n')
      fi=teta*pi/180;
      mytitle='\eta';
    elseif strcmp(ringtype,'2')
      fi=hta2*pi/180;
      mytitle='hta2';
    elseif strcmp(ringtype,'3')
      fi=hta3*pi/180;
      mytitle='hta3';
    elseif strcmp(ringtype,'4')
      fi=hta4*pi/180;
      mytitle='hta4';
    else
      error('Ring type for plotting not found')
    end

    tfi=0:pi/400:2*pi;
    r2=r1+dr;

    if ind==1
      plot(r1.*cos(tfi),r1.*sin(tfi),'k'); %plot ring with inner radius
    end
    plot(r2.*cos(tfi),r2.*sin(tfi),'k'); %plot ring with outer radius

%    xy0=zeros(size(fi));   
    X=[r1.*cos(-fi+pi/2); r2.*cos(-fi+pi/2)];
    Y=[r1.*sin(-fi+pi/2); r2.*sin(-fi+pi/2)];
    plot(X,Y,'k')

    f=(0:tickdeltafi/180*pi:2*pi);
    X=[(r2-ticklen).*cos(f); (r2+ticklen).*cos(f)];
    Y=[(r2-ticklen).*sin(f); (r2+ticklen).*sin(f)];
    plot(X,Y,'k');

    text(-0.1*dr,r1+0.6*dr,mytitle,'FontSize',20);

    r1=r2;
end


labelfi=0:deltafi:360-deltafi;
Xlabel=1.15*r1*cosd(-labelfi+90)-0.12;
Ylabel=1.05*r1*sind(-labelfi+90);
%circ=repmat('\circ',size(labelfi'));
circ=repmat('^o',size(labelfi'));
text(Xlabel,Ylabel,[int2str(labelfi') circ],'FontSize',16);

indlab2 = (labelfi>=180);
labelfi=labelfi-360;
Ylabel = Ylabel+0.12*(Ylabel>=0);
Ylabel = Ylabel-0.12*(Ylabel<0);
%circ=repmat('\circ',size(labelfi(indlab2)'));
circ=repmat('^o',size(labelfi(indlab2)'));
text(Xlabel(indlab2),Ylabel(indlab2),[int2str(labelfi(indlab2)') circ],'FontSize',16);

mytitle = strrep(workdbname,'_','\_');
text(2*dr,1.1*r1,mytitle,'FontSize',30);
text(-0.7*dr,0,title,'FontSize',20);

grid off
axis equal
%axis([-1 1 -1 1])
axis off

hold off
