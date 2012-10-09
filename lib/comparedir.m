%function comparedir(dir1,dir2)
%
% Version 1.0    
% Last modified  R O Zhurakivsky 2007-03-17 
% Created        R O Zhurakivsky 2005-09-?

atomsind
dir1=[CD.xyzdir filesep 'out6/r7b/xyz'];
%dir2='nw.out/dft_fitcd';
dir2=[CD.datadir filesep 'nw.out/step'];

%compare files by name in two directories. Comparison is done till first fullstop in name. Files with save name till
% fullstop are moved to the 'same' subdirectory


movedirname='same';

movedirname1=strcat(dir1,filesep,movedirname);
movedirname2=strcat(dir2,filesep,movedirname);
if ~isdir(movedirname1)
  mkdir(movedirname1);
end
if ~isdir(movedirname2)
  mkdir(movedirname2);
end

sfiles1 = dir(strcat(dir1,filesep,'*'));
sfiles2 = dir(strcat(dir2,filesep,'*'));


numfiles2 = size(sfiles2,1);
j=1;
for i=1:numfiles2

  dlm=strfind(sfiles2(i).name,'.');
  if isempty(dlm)
    continue
  end
  fshort2{j} = sfiles2(i).name(1:(dlm-1));
  ind2(j)=i;
  if isempty(fshort2{j})
    continue
  end  
  j=j+1;
end
numfiles2 = numel(fshort2);


numfiles = size(sfiles1,1);
for i=1:numfiles

  fname = sfiles1(i);
  dlm=strfind(sfiles1(i).name,'.');
  if isempty(dlm)
    continue
  end
  fshort = sfiles1(i).name(1:(dlm-1));
  if isempty(fshort)
    continue
  end  

  for j=1:numfiles2
    if fshort==fshort2{j}
      movefile(strcat(dir1,filesep,sfiles1(i).name),movedirname1);
      movefile(strcat(dir2,filesep,sfiles2(ind2(j)).name),movedirname2);
      break
    end
  end

end
