%rot23: make visualization of vibrational modes
%visualizes all conformations in databese from highest mode to lowest
%interactively asks for conformation number and mode number, if zero entered - quits from corresponding cycle
%if for mode entered number that is higher than number of possible modes in molecule than visualizes mode with 
%frequency is most close to entered number
%
% Version 1.0    
% Last modified  R O Zhurakivsky 2008-05-19
% Created        R O Zhurakivsky 2006-?-?

clear

atomsind
format compact
flPlot=0;
plottype=0;
sortmode=0
sortmodemode=1 %0-sort modes for visualizing ascending, 1-descending

moltype=140
theory='dftV2'  %#ok
usedpackage='Gaussian'  %#ok
onlyoriginal=1;  % process db with only original conformations
%T=420; %rThd evaporation temperature 
%T=421; %rUrd evaporation temperature 
T=456; % rAdo evaporation temperature is +183C

workdbname=['r' int2str(moltype)]   %#ok
if strcmp(usedpackage,'Gaussian')
  workdbname=[workdbname '_g'];
end
if ~strcmp(theory,'dft')
  workdbname=[workdbname '_' theory];
end
if onlyoriginal
    templ='_or';
    workdbname = [workdbname templ];
end
workdbname=[CD.dbdir filesep workdbname '.mat'] %#ok


load(workdbname,'workdb')
nframes=6;
nmodes=3*workdb(1).atomnum-6; %number of modes
color='b';

recnum=numel(workdb);
for i=1:recnum

    if isfield(workdb(i).gaussian,'T')
      tind = find(workdb(i).gaussian.T==T);
      if isempty(tind), error(['Desired temperature T=' num2str(T) 'K is not found']), end  
    else %???
      tind=1;
    end

    if isfield(workdb(i).gaussian,'MP2_6311__Gdp') && isfield(workdb(i).gaussian,'GEC')
  	    energy(i) = workdb(i).gaussian.MP2_6311__Gdp + workdb(i).gaussian.GEC(tind)/CC.encoef;
    else
	    disp('needed energy fields are not found');
        sortbyenergy=0;
    end
    sdesc(i) = {workdb(i).prop.sdesc};
    tchi(i) = workdb(i).prop.tchi;
end

if sortmode==0 %without sorting
  ind=1:recnum;
elseif sortmode==1 %sort by sdesc
  [xxx,ind]=sort(sdesc);
elseif sortmode==2 %sort by energy
  [xxx,ind]=sort(energy);
elseif sortmode==3 %sort by chi
  [xxx,ind]=sort(tchi);
elseif sortmode==4 %select
  ss=['AcbaA';'BcbaA';'AccaA';'AcbaS';'BccaS'];
  ind=[];	
  for i=1:size(ss,1)
    ind(end+1)=strcmpcellar(sdesc,ss(i,:))
  end
else 
  error('invalid sort node');
end

for i=ind

  recnumber=input('Record number?: ');
  if ~recnumber
    break
  elseif ~isempty(recnumber)
    i=recnumber;
  end


  ms0=workdb(i);


%  ms0=createzmt(ms0);
%  ms0=zmt2xyz(ms0);
%keyboard
  msfreq=freqrestruct(ms0,moltype);


if 0
    msfreq=identmol(msfreq,moltype);
    order=[16 10 12]; %'pO4' 'pC1' 'pC2'
    ind(1)=find(msfreq.pind==order(1));
    ind(2)=find(msfreq.pind==order(2));
    ind(3)=find(msfreq.pind==order(3));
    Ifr=[msfreq.x(ind) msfreq.y(ind) msfreq.z(ind)]


    ind(1)=find(ms0.pind==order(1));
    ind(2)=find(ms0.pind==order(2));
    ind(3)=find(ms0.pind==order(3));
    I=[ms0.x(ind) ms0.y(ind) ms0.z(ind)]

	X0=[msfreq.x msfreq.y msfreq.z]'
    XX=I*inv(Ifr)*X0
    X=[ms0.x ms0.y ms0.z]'
end


  clf;
  hold on
  grid off
  set(gca,'DataAspectRatio',[1 1 1])
  xlabel('X');
  ylabel('Y');
  zlabel('Z');
  cla

  if plottype==1
    light;
    light('Position',[-1 -1 -2]);
    [x,y,z] = sphere(20);
  end

  I1=strcmpcellar(msfreq.labels,'H');
  I2=strcmpcellar(msfreq.labels,'O');
  I3=strcmpcellar(msfreq.labels,'C');
  I4=strcmpcellar(msfreq.labels,'N');
  hp=plot3(msfreq.x(I1),msfreq.y(I1),msfreq.z(I1),'go',...
  		msfreq.x(I2),msfreq.y(I2),msfreq.z(I2),'ro',...
  		msfreq.x(I3),msfreq.y(I3),msfreq.z(I3),'co',...
  		msfreq.x(I4),msfreq.y(I4),msfreq.z(I4),'bo');
  hl=plot3([msfreq.x(msfreq.btA)';msfreq.x(msfreq.btB)'],... 
  		[msfreq.y(msfreq.btA)';msfreq.y(msfreq.btB)'],...
  		[msfreq.z(msfreq.btA)';msfreq.z(msfreq.btB)'],[color,'-']);

  axis manual
  lim.x=[min(msfreq.x) max(msfreq.x)];
  lim.y=[min(msfreq.y) max(msfreq.y)];
  lim.z=[min(msfreq.z) max(msfreq.z)];
  range.x=lim.x(2)-lim.x(1);
  range.y=lim.y(2)-lim.y(1);
  range.z=lim.z(2)-lim.z(1);
  dscale=0.1;
  set(gca,'XLim',[lim.x(1)-dscale*range.x lim.x(2)+dscale*range.x]);
  set(gca,'YLim',[lim.y(1)-dscale*range.y lim.y(2)+dscale*range.y]);
  set(gca,'ZLim',[lim.z(1)-dscale*range.z lim.z(2)+dscale*range.z]);

  ht1=text(msfreq.x+0.1,msfreq.y+0.1,msfreq.z+0.1,num2str(msfreq.ind'));
%  ht2=text(msfreq.x+0.2,msfreq.y+0.2,msfreq.z+0.2,strcat('p',int2str(pind.ind(msfreq.pind))));

  rotate3d on


  dfi=2*pi/nframes;
  if sortmodemode==0 %ascending
     modes=[1:+1:nmodes];
  elseif sortmodemode==1 %descending
     modes=[nmodes:-1:1];
  else
     error('unknown modes sort mode specified');
  end 
  
  for mode=modes

    modenumber=input('Mode number?: ');
    if ~modenumber
      break
    elseif ~isempty(modenumber)
      if modenumber>nmodes
        [xxx,I]=min(abs(msfreq.freq-modenumber))
		mode=I;
	  else
        mode=modenumber;
	  end
    end

    title([ms0.desc '/' ms0.prop.sdesc ': mode #' int2str(mode) ', frequency ' num2str(msfreq.freq(mode),4) 'cm^-1; intensity ' num2str(msfreq.intencity(mode),4) 'KM/Mole; force constant ' num2str(msfreq.forceconst(mode)) 'mDyne/A']);

    for k = 1:nframes
      
      fi=(k-1)*dfi;
%      cla
      msfreq.x = ms0.x + msfreq.dx(mode,:)'*sin(fi);
      msfreq.y = ms0.y + msfreq.dy(mode,:)'*sin(fi);
      msfreq.z = ms0.z + msfreq.dz(mode,:)'*sin(fi);

%      h=plotmol(msfreq);
if plottype==0
      if ~isempty(I1); set(hp(1),'XData',msfreq.x(I1),'YData',msfreq.y(I1),'ZData',msfreq.z(I1)); end
      if ~isempty(I2); set(hp(2),'XData',msfreq.x(I2),'YData',msfreq.y(I2),'ZData',msfreq.z(I2)); end
      if ~isempty(I3); set(hp(3),'XData',msfreq.x(I3),'YData',msfreq.y(I3),'ZData',msfreq.z(I3)); end
      if ~isempty(I4); set(hp(4),'XData',msfreq.x(I4),'YData',msfreq.y(I4),'ZData',msfreq.z(I4)); end
	  for h=1:numel(hl)
 	    set(hl(h),'XData',[msfreq.x(msfreq.btA(h)),msfreq.x(msfreq.btB(h))],...
				  'YData',[msfreq.y(msfreq.btA(h)),msfreq.y(msfreq.btB(h))],...
				  'ZData',[msfreq.z(msfreq.btA(h)),msfreq.z(msfreq.btB(h))]);
	  end
elseif plottype==1
  for ii = 1:mol.atomnum
    switch mol.labels(ii)
    case  'H', color = [0.7 0.7 0.7]; r = 0.3;
    case  'C', color = [0.3 0.3 1.0]; r = 0.5;
    case  'O', color = [0.3 1.0 0.3]; r = 0.5;
    case  'N', color = [1.0 0.3 1.0]; r = 0.4;
    case  'F', color = [1.0 0.0 0.5]; r = 0.4;
    case  'P', color = [1.0 1.0 0.0]; r = 0.6;
    otherwise, color = [1.0 0.0 0.0]; r = 0.5;
    end
    surface('XData',mol.x(ii) + r*x,'YData',mol.y(ii) + r*y,...
        'ZData',mol.z(ii) + r*z,'FaceColor',color,...
        'EdgeColor','none','FaceLighting','gouraud')
  end

  color = [0.7 0.7 0.7];
  if isfield(mol,'btA')	
  %plot bonds using bond table if exists
      [x,y,z] = cylinder(0.1);
      for I=1:numel(mol.btA)

        v=[mol.x(mol.btB(I))-mol.x(mol.btA(I)),...
           mol.y(mol.btB(I))-mol.y(mol.btA(I)),...
           mol.z(mol.btB(I))-mol.z(mol.btA(I))];

        Oz=v/norm(v);
        vi=[1 1 1];
        Ox=cross(Oz,vi/norm(vi));
        Oy=cross(Oz,Ox);

        r=[x(1:end); y(1:end); z(1:end)*norm(v)];
        rr=[Ox;Oy;Oz]'*r;
        xx=reshape(mol.x(mol.btA(I))+rr(1,:),size(x));
        yy=reshape(mol.y(mol.btA(I))+rr(2,:),size(x));
        zz=reshape(mol.z(mol.btA(I))+rr(3,:),size(x));

        surface('XData',xx,'YData',yy,'ZData',zz,...
				'FaceColor',color,'EdgeColor','none','FaceLighting','gouraud');

     end
  end

end

      M(k)=getframe;
    end
    movie(M,1,12);

  end

end
