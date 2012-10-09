%paste molecule properties from Excel file
%
% Version 1.0    
% Last modified  R O Zhurakivsky 2005-09-20
% Created        R O Zhurakivsky 2005-09-?

global pind

clear 
format compact
%grid on
atomsind

xlsfile = 'rib96_mainval.xls';

load 'dbrib.mat'

recnum=numel(dbrib);
for i=1:recnum
    desc(i) = {dbrib(i).desc};
end

[numeric,txt]=xlsread(xlsfile,'rib96_mainval');

ldesc=txt(3:end,3);
lnwfilename=txt(3:end,1);
lnwworktitle=txt(3:end,2);

lSPenergy=numeric(1:end,1);
lGOenergy=numeric(1:end,2);
lDM=numeric(1:end,3);
lZPE=numeric(1:end,4);

for i=1:numel(ldesc)
  j = strcmpcellar(desc,ldesc(i));
  if j
    buf = lnwfilename(i);
    dbrib(j).nwfilename=buf{1};
    buf = lnwworktitle(i);
    dbrib(j).nwworktitle=buf{1};
    dbrib(j).SP.energy = lSPenergy(i);
    dbrib(j).GO.energy = lGOenergy(i);
    dbrib(j).DM = lDM(i);
    dbrib(j).ZPE = lZPE(i);
  end
end

save 'dbrib.mat' dbrib
