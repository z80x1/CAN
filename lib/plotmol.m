function hFig=plotmol(mol,color,plottype,fl_numbersout,fl_controls,hA)
% color - bonds color
% plottype : 0 - old style (lines) / 1 - new style (spheres & cylynders)
% some code is from DRAWBDP by Joe Hicklin
%
% Version 1.0    
% Last modified  R O Zhurakivsky 2007-10-07
% Created        R O Zhurakivsky 2005-09-?

global flPlot
global pind;

atomsind;
pindsdef;

if flPlot~=0
    return
end

if nargin<1
    error('no molecular structure specified');
end
if nargin<2
    color='b';
end
if nargin<3
    plottype = 1;
end
if nargin<4
    fl_numbersout = 1;
end
if nargin<5
    fl_controls = 0;
end


if nargin<6
    scrsz = get(0, 'ScreenSize'); % размеры экрана монитора
    wscr = 640; % ширина окна
    hscr = 480; % высота окна
    waxe = wscr-40; % ширина осей
    haxe = hscr-40; % высота осей
    pl = (scrsz(3)-wscr)/2; % левый край фигуры
    pu = (scrsz(4)-hscr)/2; % верхний край фигуры

    %critd = 1.6;
    %handle=figure('Color','White');
    %set(gcf,'Color','white');

    hFig = figure('Position',[pl pu wscr hscr],'NumberTitle', 'off','Name', mol.desc);
    hA = gca;
    set(hA,'Position',[0 0 0.8 1],'visible','off','DataAspectRatio',[1 1 1])
    
    %set(hFig,'MenuBar', 'none')
    %set(hFig,'WindowButtonDownFcn', 'ChangePointer(hFig,mol);')
    set(hFig,'KeyPressFcn',@rotatefig);
    set(hFig, 'Renderer', 'opengl'); %default for complex graphics

    cameratoolbar;
else
    hFig = -1;
    set(hA,'visible','off','DataAspectRatio',[1 1 1])
end

fl_drawaxes=0;

hold on
grid off
xlabel('X');
ylabel('Y');
zlabel('Z');

% Information for all buttons
labelColor=[0.8 0.8 0.8];
top=0.95;
bottom=0.05;
left=0.75;
yInitLabelPos=0.90;
left=0.75;
labelWid=0.20;
labelHt=0.05;
btnWid=0.20;
btnHt=0.05;
% Spacing between the label and the button for the same command
btnOffset=0.003;
% Spacing between the button and the next command's label
spacing=0.05;

% The close button.
if fl_controls==1
    uicontrol( ...
        'Style','pushbutton', ...
        'Units','normalized', ...
        'Position',[left bottom btnWid 2*btnHt], ...
        'String','Close', ...
        'Callback','close(gcf)');
end

if plottype==0 %lines & points
  i1=strcmpcellar(mol.labels,'H');
  i2=strcmpcellar(mol.labels,'O');
  i3=strcmpcellar(mol.labels,'C');
  i4=strcmpcellar(mol.labels,'N');
  i5=strcmpcellar(mol.labels,'P');
  i6=strcmpcellar(mol.labels,'F');
  plot3(mol.x(i1),mol.y(i1),mol.z(i1),'go',mol.x(i2),mol.y(i2),mol.z(i2),'ro',...
        mol.x(i3),mol.y(i3),mol.z(i3),'co',mol.x(i4),mol.y(i4),mol.z(i4),'bo',...
        mol.x(i5),mol.y(i5),mol.z(i5),'yo',mol.x(i6),mol.y(i6),mol.z(i6),'mo');

  if fl_numbersout
      if ~isfield(mol,'pind') 
        for I = 1:mol.atomnum
          text(mol.x(I)+0.1,mol.y(I)+0.1,mol.z(I)+0.1,num2str(mol.ind(I)));
        end
      else
    %    text(mol.x(i1)+0.1,mol.y(i1)+0.1,mol.z(i1)+0.1,num2str(mol.ind(i1)),'Color','green');
    %    text(mol.x(i2)+0.1,mol.y(i2)+0.1,mol.z(i2)+0.1,num2str(mol.ind(i2)),'Color','red');
    %    text(mol.x(i3)+0.1,mol.y(i3)+0.1,mol.z(i3)+0.1,num2str(mol.ind(i3)),'Color','cyan');
    %    text(mol.x(i4)+0.1,mol.y(i4)+0.1,mol.z(i4)+0.1,num2str(mol.ind(i4)),'Color','blue');
    %    text(mol.x(i5)+0.1,mol.y(i5)+0.1,mol.z(i5)+0.1,num2str(mol.ind(i5)),'Color','yellow');
    %    text(mol.x(i6)+0.1,mol.y(i6)+0.1,mol.z(i6)+0.1,num2str(mol.ind(i6)),'Color','magenta');
        text(mol.x+0.1,mol.y+0.1,mol.z+0.1,num2str(reshape(mol.ind,numel(mol.ind),1)),'Color','black');

        corinds=mol.pind>0;
    %    pp=repmat('p',single(mol.atomnum),1);
        pp=repmat('p',single([numel(find(mol.pind>0)) 1]));
        
        text(mol.x(corinds)+0.3,mol.y(corinds)+0.3,mol.z(corinds)+0.3,[pp,int2str(pind.ind(mol.pind(corinds))')],'Color','blue');
      end
  end

  if isfield(mol,'btA') 
  %plot bonds using bond table if exists
      for I=1:numel(mol.btA)
        plot3([mol.x(mol.btA(I)) mol.x(mol.btB(I))],... 
            [mol.y(mol.btA(I)) mol.y(mol.btB(I))],...
            [mol.z(mol.btA(I)) mol.z(mol.btB(I))],[color,'-']);
      end
  else
      for I = 1:mol.atomnum,
          for J = I+1:mol.atomnum,
              if norm([mol.x(J)-mol.x(I) mol.y(J)-mol.y(I) mol.z(J)-mol.z(I)]) < bondlen(mol.labels(I),mol.labels(J))
                  plot3([mol.x(I) mol.x(J)], [mol.y(I) mol.y(J)], [mol.z(I) mol.z(J)],[color,'-']);
              end
          end
      end
  end

elseif plottype==1 %balls & cylinders
  light;
  light('Position',[-1 -1 -2]);
  [x,y,z] = sphere(20);
  for ii = 1:mol.atomnum
    switch mol.labels{ii}
%    case  'H', color = [0.7 0.7 0.7]; r = 0.3;
%    case  'C', color = [0.3 0.3 1.0]; r = 0.5;
%    case  'O', color = [0.3 1.0 0.3]; r = 0.5;
%    case  'N', color = [1.0 0.3 1.0]; r = 0.4;
%    case  'F', color = [1.0 0.0 0.5]; r = 0.4;
%    case  'P', color = [1.0 1.0 0.0]; r = 0.6;
%    otherwise, color = [1.0 0.0 0.0]; r = 0.5;
    case  'H', color = [0.7 0.7 0.7]; r = 0.3;
    case  'C', color = 'green'; r = 0.5;
    case  'O', color = 'red'; r = 0.5;
    case  'N', color = 'blue'; r = 0.4;
    case  'F', color = 'magenta'; r = 0.4;
    case  'P', color = 'yellow'; r = 0.6;
    case  'Na', color = 'magenta'; r = 0.6;
    case  'Cl', color = [0.2 0.2 0.2]; r = 0.7;
    otherwise, color = [0.0 0.0 0.0]; r = 0.5;
    end
    fids(ii) = surface('XData',mol.x(ii) + r*x,...
        'YData',mol.y(ii) + r*y,...
        'ZData',mol.z(ii) + r*z,...
        'FaceColor',color,...
        'EdgeColor','none',...
        'FaceLighting','gouraud', 'ButtonDownFcn','disp([''Click!'']);');
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
        Ox=cross(Oz,vi);
        Ox=Ox/norm(Ox);
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

    Hblist_ind=[];
    if isfield(mol,'HBlist') && ~isempty(mol.HBlist) % plots Hbonds specified in molecular structure

    	Hblist_ind = mol.HBlist;

    elseif isfield(mol,'AIM') && ~isempty(mol.AIM) % plots only typical 3 or 4 atoms interactions
        for I=1:size(mol.AIM.pinds,1)
            bond = mol.AIM.pinds(I,:);
            if bond(2)*bond(3)
                atmFirst=find( mol.pind==mol.AIM.pinds(I,2) .* (mol.nfrag==mol.AIM.nfrag(I,2)) );
                atmSecond=find( mol.pind==mol.AIM.pinds(I,3) .* (mol.nfrag==mol.AIM.nfrag(I,3)) );
                Hblist_ind(end+1,:) = [atmFirst atmSecond];
%                Hblist_ind(end+1,:) = [find(mol.pind==mol.AIM.pinds(I,2)) find(mol.pind==mol.AIM.pinds(I,3))'];
            end
        end
	end

    R = 0.1; %cylinder radius
    N = 20;  %number of points around the circumference  
%    [x,y,z] = cylinder(R,N); %unit cylinder
	[x,y,z] = sphere(N);
    
    for I=1:size(Hblist_ind,1)

        i1=Hblist_ind(I,1);
   	    i2=Hblist_ind(I,2);
%        i1=find(strcmpcellar(pind.labels,Hblist{I,1})==mol.pind);
%        i2=find(strcmpcellar(pind.labels,Hblist{I,2})==mol.pind);
        v=[mol.x(i2)-mol.x(i1),mol.y(i2)-mol.y(i1),mol.z(i2)-mol.z(i1)]; %bond vector

        Oz=v/norm(v);
        vi=[1 1 1];
        Ox=cross(Oz,vi);
        Ox=Ox/norm(Ox);
        Oy=cross(Oz,Ox); %Ox, Oy, Oz - orts by axes

  %      r=[x(1:end); y(1:end); z(1:end)*norm(v)/4];
%        r=[x(1:end); y(1:end); 0.08 * z(1:end)*norm(v)]; %cylinder with height 1/20
        r = 0.1 * [x(1:end); y(1:end); z(1:end)]; %sphere with radius 0.1
        rr=[Ox;Oy;Oz]'*r;
        xx=reshape(rr(1,:),size(x));
        yy=reshape(rr(2,:),size(x));
        zz=reshape(rr(3,:),size(x));

        for icyl=0:0.1:0.9 %set number and step of small cylinders
            %???maybe small spheres will be better choice than cylinders
            surface('XData',mol.x(i1)+v(1)*icyl+xx,...
            		'YData',mol.y(i1)+v(2)*icyl+yy,...
            		'ZData',mol.z(i1)+v(3)*icyl+zz,...
                    'FaceColor',color, 'EdgeColor','none', 'FaceLighting','gouraud');
        end
    end

    %h=rotate3d;
    %set(h,'ButtonDownFilter',@mycallback);
    %set(h,'Enable','on');
    %set(h,'RotateStyle','orbit');
    set(fids,'ButtonDownFcn','disp(''This executes'')');
    set(fids,'Tag','DoNotIgnore');
end %plottype==1

if isfield(mol,'desc')
  title(mol.desc);
end


if fl_drawaxes %draw axes
  maxX = max(abs(mol.x)); maxX=max([1.2*maxX 1]);
  maxY = max(abs(mol.y)); maxY=max([1.2*maxY 1]);
  maxZ = max(abs(mol.z)); maxZ=max([1.2*maxZ 1]);
  plot3([0 maxX],[0 0],[0 0],'k');
  plot3([0 0],[0 maxY],[0 0],'k');
  plot3([0 0],[0 0],[0 maxZ],'k');
  text([1.05*maxX 0 0],[0 1.05*maxY 0],[0 0 1.05*maxZ],['x';'y';'z'],'Color','black');
end

%???
if hFig > -1
    set(hFig,'Position',[pl pu wscr hscr]); %zhr100706
end

return

function rotatefig(src,evnt)

[az,el] = view;
if strcmp(evnt.Key,'j') || strcmp(evnt.Key,'leftarrow')
    view(az-1,el);
elseif strcmp(evnt.Key,'k') || strcmp(evnt.Key,'rightarrow')
    view(az+1,el);
elseif strcmp(evnt.Key,'i') || strcmp(evnt.Key,'uparrow')
    view(az,el-1);
elseif strcmp(evnt.Key,'m') || strcmp(evnt.Key,'downarrow')
    view(az,el+1);
end
return 
