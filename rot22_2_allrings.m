%rot22_2: create "conformations rings" for beta, gamma, epsilon, chi and P angles at one picture
%
% Version 1.0    
% Last modified  R O Zhurakivsky 2008-07-31
% Created        R O Zhurakivsky 2006-?-?

atomsind
format compact


%-------------------------------------------
%moltypes= [9 151 14 13 15 17 19 21 270 521 11 12 16 140 240 340 540 110 152 310 352 353 390 391 161 162] %#ok
%theories =[3   2  1  3  3  3  3  3   2   2  2  2  2   2   2   2   2   3   3   3   3   3   3   3   3   3];
moltypes= [312] %#ok
theories =[3];
postfix = ['bis'];

usedpackage='Gaussian'  %#ok
%theory='dft'          %#ok
onlyoriginal=1  %#ok %process db with only original conformations

onlynorth = 0;
onlysouth = 0;
%-------------------------------------------

theories=strcat ( 'dftV', int2str(theories'));
theories=strcat ( theories, postfix);

[tbeta,tgamma,tdelta,tepsilon,tchi,tteta,teta,P,energy]=deal(zeros(0));
Pdeg=zeros(0);
[hta2,hta3,hta4]=deal(zeros(0));

if numel(moltypes)>1
    fl_multiplemol=1;
else
    fl_multiplemol=0;
end

for mind=1:numel(moltypes)
    moltype=moltypes(mind);
    theory=theories(mind,:);

    title='';
    if ~fl_multiplemol
    if moltype==7
        title='ryb'
    elseif moltype==8
        title='dryb'
    elseif moltype==11
        title='d6AC'
    elseif moltype==16
        title='dGuo'
    elseif moltype==15
        title='dAdo'
    elseif moltype==9
        title='dCyd'
    elseif moltype==12
        title='dUrd'
    elseif moltype==13
        title='dThd'
    elseif moltype==14
        title='6AC'
    elseif moltype==110
        title='d8AA'
    elseif moltype==112
        title='d4A'
    elseif moltype==140
        title='Ado'
    elseif moltype==152
        title='d8OA'
    elseif moltype==161
        title='dP'
    elseif moltype==162
        title='d8AP'
    elseif moltype==212
        title='d4C'
    elseif moltype==240
        title='Cyd'
    elseif moltype==310
        title='d8AG'
    elseif moltype==312
        title='d4G'
    elseif moltype==340
        title='Guo'
    elseif moltype==352
        title='d8OG'
    elseif moltype==353
        title='dm7G'
    elseif moltype==390
        title='dXao'
    elseif moltype==391
        title='dIno'
    elseif moltype==412
        title='d4T'
    elseif moltype==413
        title='AZT'
    elseif moltype==512
        title='d4U'
    elseif moltype==521
        title='Br5dUrd'
    elseif moltype==523
        title='Cl5dUrd'
    elseif moltype==524
        title='F5dUrd'
    elseif moltype==540
        title='Urd'
    end
    end

    if any(moltype==[8])
        ringtypes=['b'; 'g'; 'd'; 'e'; 'p'];
    elseif any(moltype==[7])
        ringtypes=['t'; 'n'; 'b'; 'g'; 'd'; 'e'; 'p'];
    elseif any(moltype==[112,212,312,412,512])
        ringtypes=['b'; 'g'; 'c'; ]; %'p';];
    elseif any(moltype==[11,15,16,9,12,13,521,161,162,110,152,310,352,353,390,391,413,523,524])
        ringtypes=['b'; 'g'; 'd'; 'e'; 'c'; 'p'];
    %    ringtypes=['b'; 'g'; 'd'; 'e'; 'c'; 'p'; '2'; '3'; '4'];
    elseif any(moltype==[140,240,340,540,14])
        ringtypes=['t'; 'n'; 'b'; 'g'; 'd'; 'e'; 'c'; 'p'];
    end

    %inclset={'AabcA'; 'EabcA'; 'JabcA'}

    workdbname=['r' int2str(moltype)];
    if strcmp(usedpackage,'Gaussian')
      workdbname=[workdbname '_g'];
    end
    if ~strcmp(theory,'dft')
      workdbname=[workdbname '_' theory];
    end
    if onlyoriginal
        workdbname = [workdbname '_or'];
    end
    workdbfname=[CD.dbdir filesep workdbname '.mat'] %#ok

    if exist(workdbfname,'file')==2
      load(workdbfname,'workdb');
      disp(['Loaded ' workdbname]);
    else
      error('rot22_2:dbnexist','Database doesn''t exists!');
    end

    %ringtype='p' % b, g, d, e, c, p, t, n
disp(['ringtypes: ' ringtypes'])

    recnum=numel(workdb);
    for i=1:recnum
      if workdb(i).new=='Y'

        ms0=workdb(i);

        if exist('inclset','var') && ~isempty(inclset)
           if isempty(intersect(ms0.prop.sdesc,inclset))
             continue
           end
        end

        if onlynorth && workdb(i).prop.Pdeg > 90 && workdb(i).prop.Pdeg < 270
            continue
        end
        if onlysouth && (workdb(i).prop.Pdeg < 90 || workdb(i).prop.Pdeg > 270)
            continue
        end

        tbeta(end+1)=ms0.prop.tbeta;
        tgamma(end+1)=ms0.prop.tgamma;

        if isfield(ms0.prop, 'tdelta')
	        tdelta(end+1)=ms0.prop.tdelta;
	    end
        if isfield(ms0.prop, 'tepsilon')
	        tepsilon(end+1)=ms0.prop.tepsilon;
	    end
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

        hta2(end+1)=ms0.prop.hta2;
        hta3(end+1)=ms0.prop.hta3;
        if isfield(ms0.prop,'hta4')
	        hta4(end+1)=ms0.prop.hta4;
	    end

      end
    end
end %for mind

cla
hold on

%ringtype=lower(ringtype);
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

    r2=r1+dr;

    tfi=0:pi/400:2*pi;
    if ind==1
      plot(r1.*cos(tfi),r1.*sin(tfi),'k'); %plot ring with inner radius
    end
    plot(r2.*cos(tfi),r2.*sin(tfi),'k'); %plot ring with outer radius

%    xy0=zeros(size(fi));   

    tfi=(0:tickdeltafi/180*pi:2*pi);
    X=[(r2-ticklen).*cos(tfi); (r2+ticklen).*cos(tfi)];
    Y=[(r2-ticklen).*sin(tfi); (r2+ticklen).*sin(tfi)];
    plot(X,Y,'k');

    text(-0.1*dr,r1+0.6*dr,mytitle,'FontSize',20);

    X=[r1.*cos(-fi+pi/2); r2.*cos(-fi+pi/2)];
    Y=[r1.*sin(-fi+pi/2); r2.*sin(-fi+pi/2)];
    plot(X,Y,'k')

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

if ~fl_multiplemol
    mytitle = strrep(workdbname,'_','\_');
    text(2*dr,1.1*r1,mytitle,'FontSize',30);
else
    mytitle = int2str(moltypes');
    text(1.3*r1,-0.0*r1,mytitle,'FontSize',15);
end
text(-0.7*dr,0,title,'FontSize',20);

grid off
axis equal
%axis([-1 1 -1 1])
axis off

hold off


disp(['beta g+ min: ' num2str(min(tbeta(tbeta>0 & tbeta<120)))]);
disp(['beta g- min: ' num2str(min(tbeta(tbeta>-120 & tbeta<0)))]);
disp(['beta t  min: ' num2str(min([tbeta(tbeta>120) tbeta(tbeta<-120)+360]))]);
disp(['beta g+ max: ' num2str(max(tbeta(tbeta>0 & tbeta<120)))]);
disp(['beta g- max: ' num2str(max(tbeta(tbeta>-120 & tbeta<0)))]);
disp(['beta t  max: ' num2str(max([tbeta(tbeta>120) tbeta(tbeta<-120)+360]))]);
disp(['beta g+ mean: ' num2str(mean(tbeta(tbeta>0 & tbeta<120)))]);
disp(['beta g- mean: ' num2str(mean(tbeta(tbeta>-120 & tbeta<0)))]);
disp(['beta t  mean: ' num2str(mean([tbeta(tbeta>120) tbeta(tbeta<-120)+360]))]);
disp(['beta g+ std: ' num2str(std(tbeta(tbeta>0 & tbeta<120)))]);
disp(['beta g- std: ' num2str(std(tbeta(tbeta>-120 & tbeta<0)))]);
disp(['beta t  std: ' num2str(std([tbeta(tbeta>120) tbeta(tbeta<-120)+360]))]);

disp(['gamma g+ mean: ' num2str(mean(tgamma(tgamma>0 & tgamma<120)))]);
disp(['gamma g- mean: ' num2str(mean(tgamma(tgamma>-120 & tgamma<0)))]);
disp(['gamma t  mean: ' num2str(mean([tgamma(tgamma>120) tgamma(tgamma<-120)+360]))]);
disp(['gamma g+ std: ' num2str(std(tgamma(tgamma>0 & tgamma<120)))]);
disp(['gamma g- std: ' num2str(std(tgamma(tgamma>-120 & tgamma<0)))]);
disp(['gamma t  std: ' num2str(std([tgamma(tgamma>120) tgamma(tgamma<-120)+360]))]);

disp(['epsilon g+ mean: ' num2str(mean(tepsilon(tepsilon>0 & tepsilon<120)))]);
disp(['epsilon g- mean: ' num2str(mean(tepsilon(tepsilon>-120 & tepsilon<0)))]);
disp(['epsilon t  mean: ' num2str(mean([tepsilon(tepsilon>120) tepsilon(tepsilon<-120)+360]))]);
disp(['epsilon g+ std: ' num2str(std(tepsilon(tepsilon>0 & tepsilon<120)))]);
disp(['epsilon g- std: ' num2str(std(tepsilon(tepsilon>-120 & tepsilon<0)))]);
disp(['epsilon t  std: ' num2str(std([tepsilon(tepsilon>120) tepsilon(tepsilon<-120)+360]))]);




