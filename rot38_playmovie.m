%rot38_playmovie.m:  
%
% Version 1.0    
% Last modified  R O Zhurakivsky 2010-07-06
% Created        R O Zhurakivsky 2010-07-06

tic
clear 
format compact

atomsind

load('D:\_diplom\CAN\db\r13_AabcA_g_dftV3.mat')
mov = avifile('example6.avi');

recnum=numel(workdb);
for i=1:recnum
    param(i) = workdb(i).param;
end

param(param<0)=param(param<0)+360;
[xxx,ind]=sort(param);


%try
clear m;
step=0;
    for i=ind
step=step+1;
    	mol = workdb(i);

        disp([int2str(step) ' ' num2str(mol.param)]);

    	h = plotmol(mol,'b',1,1,0);
        curaxis = axis;
        x = curaxis(1) + (curaxis(2)-curaxis(1))/20;
        y = curaxis(4) - (curaxis(4)-curaxis(3))/20;
        h_t = text(x,y,['\chi = ' num2str(mol.param,3) ' град']);
        set(h_t,'FontSize',14,'FontWeight','bold');

    	m(step) = getframe(h,[0 0 640 480]);

%if i==32, keyboard, end
       	mov = addframe(mov,m(step));

    	close(h);

    end

    mov=close(mov);

%    movie(m,-10);


%catch
%    warning('rot38_playmovie:exception occured');
%	mov=close(mov);
%end

