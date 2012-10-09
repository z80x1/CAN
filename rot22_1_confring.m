%rot22: create "conformations rings" for beta, gamma, epsilon, chi and P angles
%
% Version 1.0    
% Last modified  R O Zhurakivsky 2008-07-31
% Created        R O Zhurakivsky 2006-?-?

atomsind
format compact

moltype=112
workdbname=[CD.dbdir filesep 'r' int2str(moltype) '_g_dftV3bis_or.mat']
%workdbname=[CD.dbdir filesep 'r8_g_or.mat']
ringtype='e' % b, g, d, e, c, p, t, n
onlynorth = 0
onlysouth = 0

load(workdbname,'workdb')

recnum=numel(workdb);

[tbeta,tgamma,tdelta,tepsilon,tchi,tteta,teta,P,energy]=deal(zeros(0));
Pdeg=zeros(0);
for i=1:recnum
  if workdb(i).new=='Y'

    ms0=workdb(i);

    if onlynorth && workdb(i).prop.Pdeg > 90 && workdb(i).prop.Pdeg < 270
        continue
    end
    if onlysouth && (workdb(i).prop.Pdeg < 90 || workdb(i).prop.Pdeg > 270)
        continue
    end

    tbeta(end+1)=ms0.prop.tbeta;
    tgamma(end+1)=ms0.prop.tgamma;
    tdelta(end+1)=ms0.prop.tdelta;
    tepsilon(end+1)=ms0.prop.tepsilon;
    if isfield(ms0.prop,'tchi')
        tchi(end+1)=ms0.prop.tchi;
    end
    if isfield(ms0.prop,'tteta')
        tteta(end+1)=ms0.prop.tteta;
    end
    if isfield(ms0.prop,'teta')
        teta(end+1)=ms0.prop.teta;
    end
%    P(end+1)=ms0.prop.P;
    Pdeg(end+1)=ms0.prop.Pdeg;
    energy(end+1)=ms0.GO.energy;
%   energy(end+1)=ms0.gaussian.MP2_6311__G2dfpd;

  end
end

denergy=(max(energy)-energy)*CC.encoef; %energy delta

ringtype=lower(ringtype);
if strcmp(ringtype,'b')
  fi=tbeta*pi/180;
elseif strcmp(ringtype,'g')
  fi=tgamma*pi/180;
elseif strcmp(ringtype,'d')
  fi=tdelta*pi/180;
elseif strcmp(ringtype,'e')
  fi=tepsilon*pi/180;
elseif strcmp(ringtype,'c')
  fi=tchi*pi/180;
elseif strcmp(ringtype,'p')
  fi=Pdeg*pi/180;
  %fi=P;
elseif strcmp(ringtype,'t')
  fi=tteta*pi/180;
elseif strcmp(ringtype,'n')
  fi=teta*pi/180;
else
  error('Ring type for plotting not found')
end

r=denergy;
maxr=max(r);

r=repmat(maxr,size(r));

x=r./sqrt(1+tan(fi).^2).*sign(fi);
y=x.*tan(fi);

x2=max(r)./sqrt(1+tan(fi).^2).*sign(fi);
y2=x2.*tan(fi);

%h=polar([zeros(size(fi));fi],[zeros(size(r));r],'k');

tfi=0:pi/400:2*pi;
tr=repmat(1,size(tfi));

cla
hold on
plot(tr.*cos(tfi),tr.*sin(tfi),'k');

xy0=zeros(size(fi));	
X=[xy0; r.*cos(-fi+pi/2)]/maxr;
Y=[xy0; r.*sin(-fi+pi/2)]/maxr;

plot(X,Y,'k')
deltafi=30; %interval of fi labeling
tickdeltafi=15; %interval of tick placing
ticklen=0.02;

labelfi=0:deltafi:360-deltafi;
Xlabel=1.15*cosd(-labelfi+90)-0.06;
Ylabel=1.05*sind(-labelfi+90);
text(Xlabel,Ylabel,int2str(labelfi'),'FontSize',16);

indlab2 = (labelfi>=180);
labelfi=labelfi-360;
Ylabel = Ylabel+0.1*(Ylabel>=0);
Ylabel = Ylabel-0.1*(Ylabel<0);
text(Xlabel(indlab2),Ylabel(indlab2),int2str(labelfi(indlab2)'),'FontSize',16);
	

fi=(0:tickdeltafi/180*pi:2*pi);
X=[(1-ticklen).*cos(fi); (1+ticklen).*cos(fi)];
Y=[(1-ticklen).*sin(fi); (1+ticklen).*sin(fi)];
plot(X,Y,'k');


grid off
axis equal
%axis([-1 1 -1 1])
axis off

switch ringtype
case 'b'
  mytitle='\beta';
  drawarc([pi/6,pi/2],'g^+');
  drawarc([5*pi/6,7*pi/6],'t');
  drawarc([-pi/6,-pi/2],'g^-');
case 'g'
  mytitle='\gamma';
  drawarc([pi/6,pi/2],'g^+');
  drawarc([5*pi/6,7*pi/6],'t');
  drawarc([-pi/6,-pi/2],'g^-');
case 'd'
  mytitle='\delta';
  drawarc([pi/6,pi/2],'g^+');
  drawarc([5*pi/6,7*pi/6],'t');
  drawarc([-pi/6,-pi/2],'g^-');
case 'e'
  mytitle='\epsilon';
  drawarc([pi/6,pi/2],'g^+');
  drawarc([5*pi/6,7*pi/6],'t');
  drawarc([-pi/6,-pi/2],'g^-');
case 't'
  mytitle='\theta';
  drawarc([pi/6,pi/2],'g^+');
  drawarc([5*pi/6,7*pi/6],'t');
  drawarc([-pi/6,-pi/2],'g^-');
case 'n'
  mytitle='\eta';
  drawarc([pi/6,pi/2],'g^+');
  drawarc([5*pi/6,7*pi/6],'t');
  drawarc([-pi/6,-pi/2],'g^-');
case 'c'
  mytitle='\chi';
  drawarc([-pi/2,pi/2],'syn');
  drawarc([pi/2,3/2*pi],'anti');
case 'p'
  mytitle='P';
  drawarc([-1/4*pi, 1/4*pi],'N');
  drawarc([ 3/4*pi, 5/4*pi],'S');
  drawarc([-3/4*pi,-1/4*pi],'W');
  drawarc([ 1/4*pi, 3/4*pi],'E');
end
%mytitle=['Gmax-G(' mytitle ')'];

if onlynorth 
    mytitle = [mytitle '_N'];
end
if onlysouth 
    mytitle = [mytitle '_S'];
end

text(-0.20,1.20,mytitle,'FontSize',60);

hold off
