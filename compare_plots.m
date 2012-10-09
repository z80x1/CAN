%load data from txt files and plot them in one figure
% 2010-0112 used for E(chi) dependencies for nucleosides


%--------- data section -------------------------------
indir = 'D:\_diplom\CAN2\data\rx21' %#ok
ext = 'res2' %#ok
%--------- data section -------------------------------

pcolor=[{'r'} {'g'} {'b'} {'k'} {'m'} {'c'} ...
        {[.8627 .0784 .2353]} {[.5412 .1686 .8863]}  ...
        {[0 .5 0]} {[.5 0 0]} {[0 0 .5]} {[.5 .5 0]} {[.5 0 .5]} {[0 .5 .5]} {[.5 .5 .5]} ...
        {[.5 .1 .9]} {[.8 .2 .2]} {[.8 .8 .2]}...
        {[.9 .4 .9]} {[.2 .4 .6]} {[.6 .4 .6]} {[.6 .2 .2]} {[.8 .2 .8]} ...
        {[.2 .8 .8]} {[.2 .8 .2]} {[.2 .2 .8]} {[.4 .9 .1]} {[.1 .3 .6]} {'y'} {[.75 .75 .75]} {[.2745 .5098 .7059]}];


sfiles = dir(strcat(indir,filesep,['*.' ext]));
numfiles = size(sfiles,1);
if ~numfiles, error('No files found!'), end

figure
hold on

h=[]; %plot handles array
legs=[]; %legends array

for f_ind=1:numfiles

    fname=sfiles(f_ind).name;
    dlm=strfind(fname,'.');
    fnameshort = fname(1:(dlm-1));

    fid=fopen([indir filesep fname],'r'); % don't use text mode - this is dangerous for Unix files ;)
    C = textscan(fid,'%f%f');
    fclose(fid);

    x = C{:,1};
    x(x<0) = x(x<0)+360;
    [x,I]= sort(x);
    y = C{:,2};
    y = y(I);

    h(end+1) = plot(x,y,['.-'],'Color', pcolor{f_ind},'MarkerSize', 6);
    legs{end+1}= fnameshort;

end

hl=legend(h, legs,'Location','NEO');
xlabel('\chi');
ylabel('energy, kcal/mol');
oldaxis=axis;
axis([0 360 0 oldaxis(4)]);
grid on
grid minor
set(gca,'XTick',[0:45:360]);
hold off
   

%A=polyfit(C{:,1},C{:,2},2)
%plot(C{:,1},A(1).*C{:,1}.^2+A(2).*C{:,1}+A(3),'.-g')
