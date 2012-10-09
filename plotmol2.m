function plotmol2(flag)
% some code is from DRAWBDP by Joe Hicklin
%!!do not save images in image window  maximize mode
%
% Version 1.0    
% Last modified  R O Zhurakivsky 2008-07-31
% Created        R O Zhurakivsky 2006-?-?

global flPlot
global pind

global sdesc ind workdb workname

pindsdef
atomsind

if flPlot~=0
    return
end

%bgcolor='w';
fl_showtext=0;
fl_plotgray=0;
fl_labelatomsbytype=0;
fl_plotlines=0;
fl_showhbonds=1;
fl_showcharges=0;

k_bondthick = 0.08;

persistent lst; % the list box

if(nargin == 0)
%starting branch
%    handle=figure('Color',bgcolor);
    set(gca,'Position',[0 0 0.8 1],'visible','off','DataAspectRatio',[1 1 1])
    cameratoolbar;
    sdesc={};

  %---------------
    workname='r391'  %#ok
    onlyoriginal=1  %#ok % process db with only original conformations
    theory='dftV3'  %#ok
  %---------------

    workdbname=[CD.dbdir filesep workname '_g'];
    if ~strcmp(theory,'dft')
      workdbname=[workdbname '_' theory];
    end
    if onlyoriginal
  	workdbname=[workdbname '_or'];
    end
    workdbname=[workdbname '.mat'] %#ok

    load(workdbname,'workdb')
    recnum=numel(workdb);
    for i=1:recnum
      sdesc(i) = {workdb(i).prop.sdesc};
    end
    [sdesc_,ind]=sort(sdesc);
    sdesc_=[sdesc_,{'atoms'}];

    lst = uicontrol('Units','normalized', ...
          'Position',[.8, .05,.19,.9],...
          'String',sdesc_,'Style','listbox','Callback','plotmol2(1)');

    if ~isfield(workdb(1).gaussian,'mcharge')
        fl_showcharges=0;
    end


else %plotting specified structure

    cla
    hold on
    grid off
    light;
    light('Position',[-1 -1 -2]);

    nm = get(lst,'String');

    val=get(lst,'Value');
  %keyboard
%    plotatoms=0;
    if val==numel(get(lst,'String'))
%      plotatoms=1;
        mol.atomnum=4;
        mol.labels=[{'H'} {'O'} {'C'} {'N'}];
        mol.x=[0;0;0;0];
      	mol.y=[0;-1;-2.2;-3.2];
        mol.z=[0;0;0;0];
        axis([-5 5 -5 5 -5 5]);
        mol.gaussian.mcharge=[0;0;0;0];
    else
        mol=workdb(ind(val));
        axis tight
    end

    [x,y,z] = sphere(20);

    for ii = 1:mol.atomnum
        if ~fl_plotgray
            switch mol.labels{ii}
            case  'H', color = [0.7 0.7 0.7]; r = 0.25;
  %          case  'C', color = [0.3 0.3 1.0]; r = 0.5;
  %          case  'O', color = [0.3 1.0 0.3]; r = 0.5;
  %          case  'N', color = [1.0 0.3 1.0]; r = 0.4;
            case  'C', color = [0 1 0]; r = 0.35;
            case  'O', color = [1 0 0]; r = 0.35;
            case  'N', color = [0 0 1	]; r = 0.3;
            otherwise, color = [1.0 0.0 0.0]; r = 0.35;
            end
        else
            switch mol.labels(ii)
            case  'H', color = [0.7 0.7 0.7]; r = 0.25;
            case  'C', color = [0.5 0.5 0.5]; r = 0.35;
            case  'O', color = [0.4 0.4 0.4]; r = 0.35;
            case  'N', color = [1.0 1.0 1.0]; r = 0.3;
            otherwise, color = [1.0 1.0 1.0]; r = 0.35;
            end
        end
      	if fl_plotlines
          	r=r*0.1;
        end

        surface('XData',mol.x(ii) + r*x,'YData',mol.y(ii) + r*y,...
                'ZData',mol.z(ii) + r*z,'FaceColor',color,...
                'EdgeColor','none','FaceLighting','gouraud')
        if fl_showtext
            text(mol.x(ii)+0.9*r,mol.y(ii)+0.9*r,mol.z(ii)+0.9*r,num2str(mol.ind(ii)),'Color',color);
            text(mol.x(ii)-0.9*r,mol.y(ii)+0.9*r,mol.z(ii)+0.9*r,['p',int2str(pind.ind(mol.pind(ii)))],'Color','white');
        end
        if fl_labelatomsbytype
  %     if mol.labels(ii)~='H'
            text(mol.x(ii)+1.3*r,mol.y(ii)+1.3*r,mol.z(ii)+1.3*r,mol.labels(ii),'Color','Black','FontSize',16);
  %     end
        end
        if fl_showcharges 
            text(mol.x(ii)+0.9*r,mol.y(ii)+0.9*r,mol.z(ii)+0.9*r,num2str(mol.gaussian.mcharge(ii),3),'Color','Black','FontSize',16);
        end
    
    end

    color = [0.7 0.7 0.7];
    if isfield(mol,'btA')	
        %plot bonds using bond table if exists
        if ~fl_plotlines

          [x,y,z] = cylinder(k_bondthick);
          for I=1:numel(mol.btA)

      	    i1=mol.btA(I); i2=mol.btB(I);
            v=[mol.x(i2)-mol.x(i1),mol.y(i2)-mol.y(i1),mol.z(i2)-mol.z(i1)];

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
        else
          for I=1:numel(mol.btA)
            plot3([mol.x(mol.btA(I)) mol.x(mol.btB(I))],... 
    			    [mol.y(mol.btA(I)) mol.y(mol.btB(I))],...
    			    [mol.z(mol.btA(I)) mol.z(mol.btB(I))],['g','-']);
          end
        end
    end

    Hblist={};
    if fl_showhbonds
        if isfield(mol,'AIM') 
        %extract information about Hbonds automatically
            for i=1:size(mol.AIM.pinds,1)
                Hblist(i,1:2)=[pind.labels(mol.AIM.pinds(i,2)) pind.labels(mol.AIM.pinds(i,3))];
            end
        else
            if strcmp(workname,'r11')
              switch nm{get(lst,'Value')}
              case 'AbcaA'
                Hblist=[{'pH31'},{'bN6'};{'pH51'},{'bN6'};{'pH12'},{'bO2'}]; %r11:AbcaA
              case 'AabcA'
                Hblist=[{'pH31'},{'bN6'};{'pH12'},{'bO2'}]; %r11:AabcA %46
              case 'EaaaS'
                Hblist=[{'pH21'},{'bO2'};{'pH21'},{'pO5'};{'pH53'},{'bO2'}]; %r11:EaaaS
              case 'EaacS'
                Hblist=[{'pH21'},{'bO2'};{'pH21'},{'pO5'};{'pH53'},{'bO2'}]; %r11:EaacS
              case 'EaacA'
                Hblist=[{'pH12'},{'bO2'};{'pH21'},{'pO5'};{'pH53'},{'bN6'}]; %r11:EaacA
              case 'FaaaA'
                Hblist=[{'pH12'},{'bO2'};{'pH21'},{'pO5'};{'pH53'},{'bN6'}]; %r11:FaaaA
              case 'GacbA'
                Hblist=[{'pH21'},{'pO5'}]; %r11:GacbA
              otherwise
                Hblist=[];
              end
            elseif strcmp(workname,'r7')
              switch nm{get(lst,'Value')}
              case 'Aaaca'
                Hblist=[{'pH22'},{'pO3'}]; %r7: #2
              case 'Acbaa'
                Hblist=[{'pH22'},{'pO3'};{'pH32'},{'pO5'}]; %r7: #2
              case 'Eabca'
                Hblist=[{'pH32'},{'pO2'};{'pH21'},{'pO5'}]; %r7: #22
              case 'Aabca'
                Hblist=[{'pH22'},{'pO3'}]; %r7: #27
              case 'Dabcb'
                Hblist=[]; %r7: #32
              case 'Eabcc'
                Hblist=[{'pH32'},{'pO2'};{'pH21'},{'pO5'}]; %r7: #24
              otherwise
                Hblist=[];
              end
            elseif strcmp(workname,'r8')
              switch nm{get(lst,'Value')}
              case {'Acca', 'Acba'}
                Hblist=[{'pH32'},{'pO5'}];
              case {'Acab', 'Bcac'}
                Hblist=[{'pH53'},{'pO3'}];
              case {'Eabc', 'Eaba', 'Gabb', 'Eabb', 'Eacc', 'Eaca', 'Gacb', 'Eacb'}
                Hblist=[{'pH21'},{'pO5'}];
              otherwise
                Hblist=[];
              end
            elseif strcmp(workname,'r9') %---!!didn't tested by AIM
              switch nm{get(lst,'Value')}
              case {'AabcA', 'AabcB'}
                Hblist=[{'bH6'},{'pO5'};{'bH6'},{'pO4'};{'bO2'},{'pH12'}];
              case {'EabcA', 'EabcC'}
                Hblist=[{'bH6'},{'pO5'};{'bH6'},{'pO4'};{'bO2'},{'pH12'}];
              case {'EaaaS', 'EaaaV'}
                Hblist=[{'pH53'},{'bO2'}];
              case {'EaabS', 'EaabV'}
                Hblist=[{'pH53'},{'bO2'}];
              case {'EaacS', 'EaacV'}
                Hblist=[{'pH53'},{'bO2'}];
              otherwise
                Hblist=[];
              end
            elseif strcmp(workname,'r12')
              switch nm{get(lst,'Value')}
              case 'EaacS'
                Hblist=[{'pH21'},{'bO2'};{'pH53'},{'bO2'};{'pH21'},{'pO5'}];
              case 'EabcA'
                Hblist=[{'pH12'},{'bO2'};{'bH6'},{'pO5'};{'pH21'},{'pO5'}];
              case 'AabcA'
                Hblist=[{'bH6'},{'pO4'};{'bH6'},{'pO5'}];
              case 'EcbcA'
                Hblist=[{'bH6'},{'pO4'}];
              case 'AcbcA'
                Hblist=[{'bH6'},{'pO4'}];
              case 'FbbcA'
                Hblist=[{'bH6'},{'pO4'};{'bH6'},{'pO5'}];
              case 'AbbcA'
                Hblist=[{'bH6'},{'pO4'}];
              case 'AbccS'
                Hblist=[{'pH31'},{'bO2'};{'pH51'},{'bO2'}];
              case 'AbbcS'
                Hblist=[{'pH31'},{'bO2'}];
              otherwise
                Hblist=[];
              end
            end
        end
    end

    [x,y,z] = cylinder(k_bondthick);
    for I=1:size(Hblist,1)

      	i1=find(strcmpcellar(pind.labels,Hblist{I,1})==mol.pind);
        i2=find(strcmpcellar(pind.labels,Hblist{I,2})==mol.pind);
        v=[mol.x(i2)-mol.x(i1),mol.y(i2)-mol.y(i1),mol.z(i2)-mol.z(i1)];

        Oz=v/norm(v);
        vi=[1 1 1];
        Ox=cross(Oz,vi);
        Ox=Ox/norm(Ox);
        Oy=cross(Oz,Ox);

  %      r=[x(1:end); y(1:end); z(1:end)*norm(v)/4];
        r=[x(1:end); y(1:end); z(1:end)*norm(v)/20];
        rr=[Ox;Oy;Oz]'*r;
        xx=reshape(rr(1,:),size(x));
        yy=reshape(rr(2,:),size(x));
        zz=reshape(rr(3,:),size(x));

        for icyl=0:2:18
            surface('XData',mol.x(i1)+v(1)*icyl/20+xx,'YData',mol.y(i1)+v(2)*icyl/20+yy,'ZData',mol.z(i1)+v(3)*icyl/20+zz,...
                    'FaceColor',color,'EdgeColor','none','FaceLighting','gouraud');
        end
    end

    if isfield(mol,'desc')
      title(mol.desc);
    end
    rotate3d on
    colormap gray

%    set(gcf,'KeyPressFcn',@rotatefig);
end

end

function rotatefig(src,evnt)

[az,el] = view;
if evnt.Key == 'j'
	view(az-1,el);
elseif evnt.Key == 'k'
	view(az+1,el);
elseif evnt.Key == 'i'
	view(az,el-1);
elseif evnt.Key == 'm'
	view(az,el+1);
end

text(0, 0, ['az=' int2str(az) ', el=' int2str(el)]);

end
