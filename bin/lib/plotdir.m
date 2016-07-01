function plotdir(indir,moltype)
%plot molecules at 'indir' in cycle
%
% Version 1.0    
% Last modified  R O Zhurakivsky 2007-03-17
% Created        R O Zhurakivsky 2006-?-?

if nargin<1
  indir='..\..\work\_yevgen\060530\'
end
if nargin<2
  moltype=0;
end

format compact
%grid on

global flPlot;
flPlot=0;

atomsind
maxind = 0;

sfiles = dir(strcat(indir,filesep,'*.xyz'));
numfiles = size(sfiles,1);
for i=1:numfiles

    delete(gca);
    dlm=strfind(sfiles(i).name,'.');
    fname = sfiles(i).name(1:(dlm-1));

    ffname = strcat(indir,filesep,sfiles(i).name); 

    ms0=loadmolxyz(ffname);
    ms0 = createbondtable(ms0);
    
    if moltype
      ms0 = identmol(ms0,moltype);
    end
    strcat('Loading file: ',ffname)

    plotmol(ms0)

%    'Press any key.'
keyboard
%    pause
end
