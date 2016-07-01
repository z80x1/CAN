%rot28_1_importAIM: import AIM data to database
%
% Version 1.0    
% Last modified  R O Zhurakivsky 2008-05-19
% Created        R O Zhurakivsky 2006-?-?

%2006-0903 - changed bond name parsing algorithm

tic
clear 
format compact

global pind
atomsind

%---------------------------------
moltype=161
usedpackage='Gaussian'
theory='dftV3'
onlyoriginal=1;  % process db with only original conformations
data='081124';
flwritefile=1 %save db file
%---------------------------------


workdbname=['r' int2str(moltype)]
ifilename = workdbname;

if usedpackage=='Gaussian'
  workdbname=[workdbname '_g'];
end
if ~strcmp(theory,'dft')
  workdbname=[workdbname '_' theory];
%  ifilename=[ifilename '_' theory];
end
ifilename=[CD.AIM filesep ifilename '_AIM.' data '.res']

if onlyoriginal
    templ='_or';
    workdbname = [workdbname templ];
end
workdbname=[CD.dbdir filesep workdbname '.mat']


fid=fopen(ifilename,'r'); % don't use text mode - this is dangerous for Unix files ;)
if fid==-1
 error(['Can''t open file ' ifilename])
end
%A=textscan(fid,'%s%d%s%s%s%s');
A=textscan(fid,'%d16%d16%s%u8%u8%s%f%f%s%f%f%s%f%f%s%f%f%s%f%f','delimiter','\t');
fclose(fid);
iConfSdesc=3; %index of field with sdesc data
iHBondsNum=4; %index of field with number of Hbonds
i1Hbondfield=6; %index of first field with Hbond description

load(workdbname,'workdb')

recnum=numel(workdb);

sdesc={};
[aind1,aind2,aind3]=deal(uint16(0));

for i=1:recnum
   sdesc{i} = workdb(i).prop.sdesc;
end

bondstotal = 0;


for i=1:numel(A{1}) %cycle though lines of input file
  ind = strcmpcellar(sdesc,A{iConfSdesc}(i));
  if isempty(ind)
    error(['Conformer ' A{iConfSdesc}(i) ' is not found between working molecule structure conformers' ]); 
  end
  
  AIM.desc={};
  AIM.ro=[];
  AIM.deltaroAIM2000=[];
  AIM.pinds=[];

  j=0;
%  skippedfields=0;
  while 1 %cycle though bonds
    j=j+1;
%    if j > A{iHBondsNum}(i)+skippedfields
    if j > 5 %!! More then 5 Hbonds are considered. !!
      break;
    end;

    bondstr=A{(j-1)*3+i1Hbondfield}{i};
    if isempty(bondstr)
%	skippedfields=skippedfields+1
	continue;
    end

    AIM.desc(j)={bondstr};
    AIM.ro(j)=A{(j-1)*3+i1Hbondfield+1}(i);
    AIM.deltaroAIM2000(j)=A{(j-1)*3+i1Hbondfield+2}(i);

    charinds = find(isstrprop(bondstr,'alpha'));
    numinds = find(isstrprop(bondstr,'digit'));
    pointind=strfind(bondstr,'.');  %instead of ... string export from Excel may contain : one

    charinds2 = charinds(charinds>numinds(1)); % charinds2(1) is first symbol of atom2 name

    atom1 = bondstr(1:charinds2(1)-1);

    atom2=bondstr(charinds2(1):pointind(1)-1);

    atom3 = bondstr(pointind(end)+1:end);

    quotesind = strfind(atom1,'''');
    if ~isempty(quotesind)
      atom1=['p' atom1(1:quotesind-1)];
    else
      atom1=['b' atom1];
    end
    aind1=strcmpcellar(pind.labels,atom1);

    quotesind = strfind(atom2,'''');
    if ~isempty(quotesind)
      atom2=[atom2(1:quotesind-1)];
    end
    diginds = find(isstrprop(atom2,'digit'));
    if isempty(diginds)
      diginds = find(isstrprop(atom1,'digit'));
      atom2=[atom2 atom1(diginds)];
    end
    atom2=[atom1(1) atom2];
    if strcmp(atom2,'pH5')
      atom2 = 'pH53';
    end
    aind2=strcmpcellar(pind.labels,atom2);

    quotesind = strfind(atom3,'''');
    if ~isempty(quotesind)
      atom3=['p' atom3(1:quotesind-1)];
    else
      atom3=['b' atom3];
    end
    aind3=strcmpcellar(pind.labels,atom3);
%    [atom1 '|' atom2 '|' atom3]

    AIM.pinds(j,1:3)=[aind1 aind2 aind3];

  end
  [{[int2str(i) ':']} AIM.desc]  %#ok
  bondstotal = bondstotal + numel(AIM.desc);
  workdb(ind).AIM = AIM;

end

disp([' Total ' int2str(bondstotal) ' of H-bonds are imported.']);

if flwritefile
  dlm=strfind(workdbname,'.');
  workdbnameold=[workdbname(1:dlm(end)-1) '~' workdbname(dlm(end):end)];
  copyfile(workdbname,workdbnameold);

  save(workdbname,'workdb')
end


toc
