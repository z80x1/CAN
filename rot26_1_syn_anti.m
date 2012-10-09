%compare geometry parameters of two sets of conformations
%
% Version 1.0    
% Last modified  R O Zhurakivsky 2008-05-19
% Created        R O Zhurakivsky 2006-?-?

clear T
clear

atomsind
format compact

workname='r11_g_dftV3_or'
workdbname=[CD.dbdir '/' workname '.mat']
xlsfile = [CD.xlsdir '/' workname '_GC' '.xls'];

infile1=[CD.dbdir '/' 'r11_anti.set']
infile2=[CD.dbdir '/' 'r11_syn.set']

fid=fopen(infile1);
i=1;
while 1 
  tline=fgetl(fid);
  if ~ischar(tline),break, end
  set1{i}=tline;
  i=i+1;
end
fclose(fid);

fid=fopen(infile2);
i=1;
while 1 
  tline=fgetl(fid);
  if ~ischar(tline),break, end
  set2{i}=tline;
  i=i+1;
end
fclose(fid);

if numel(set1)~=numel(set2)
  error('number of conformations in sets must be equal')
end
N=numel(set1);

load(workdbname,'workdb')
for i=1:numel(workdb)
  sdesc(i) = {workdb(i).prop.sdesc};
end

T(1,1)={'conf1'};
T(1,2)={'conf2'};

flist=[{'tbeta'} {'tgamma'} {'tdelta'} {'tepsilon'} {'Pdeg'} {'numax'} {'tau01'} {'tau02'} {'tau03'} {'maxtorinbase'}];
for i=1:numel(flist)
 T(1,3+3*(i-1))={[flist{i} '1']};
 T(1,4+3*(i-1))={[flist{i} '2']};
 T(1,5+3*(i-1))={['abs d' flist{i}]};
end

for i=1:N

  j1=strcmpcellar(sdesc,set1{i});
  j2=strcmpcellar(sdesc,set2{i});

  T(i+1,1)=set1(i);
  T(i+1,2)=set2(i);

  for ii=1:numel(flist)
    val1=getfield(workdb(j1).prop,flist{ii});
    val2=getfield(workdb(j2).prop,flist{ii});
    dval=val1-val2;
	dval=dval-360*fix(dval/180);

    T(i+1,3+3*(ii-1))={abs(val1)};
    T(i+1,4+3*(ii-1))={abs(val2)};
    T(i+1,5+3*(ii-1))={abs(dval)};

  end

end


xlswrite(xlsfile,T,'geom comparison','A1');
disp('Done!')
