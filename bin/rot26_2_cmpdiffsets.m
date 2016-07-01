%compare geometry parameters of two sets of conformations
%
% Version 1.0    
% Last modified  R O Zhurakivsky 2008-05-19
% Created        R O Zhurakivsky 2006-?-?

clear

atomsind
format compact

%----------------------------------------
workname='r15_g_dftV3_or'
workdbname=[CD.dbdir filesep workname '.mat']

workname2='r161_g_dftV3_or'
workdbname2=[CD.dbdir filesep workname2 '.mat']

xlsfile = [CD.xlsdir filesep workname '_r15vsr161' '.xls'];
%----------------------------------------


load(workdbname,'workdb')
lsdesc = numel(workdb(1).prop.sdesc);

workdb1 = workdb;
for i=1:numel(workdb1)
  sdesc(i) = {workdb1(i).prop.sdesc};
  GOenergy1(i) = workdb1(i).GO.energy;
  DM1(i) = workdb1(i).DM;
end
minGOenergy1 = min(GOenergy1);
dGOenergy1 = (GOenergy1-minGOenergy1)*CC.encoef;

load(workdbname2,'workdb')
workdb2 = workdb;
N2 = numel(workdb2);
for i=1:N2
  sdesc2(i) = {workdb2(i).prop.sdesc};
  GOenergy2(i) = workdb2(i).GO.energy;
  DM2(i) = workdb2(i).DM;
end
minGOenergy2 = min(GOenergy2);
dGOenergy2 = (GOenergy2-minGOenergy2)*CC.encoef;

sdesc_com = intersect(sdesc,sdesc2);

sdesc_new1 = sdesc_com;
sdesc_new2 = sdesc_com;

sdesc_np = setdiff(sdesc2,sdesc_com);
for i=1:numel(sdesc_np)
  [iprev,inext,sdesc_prev,sdesc_next] = checksimdesc(sdesc_np{i},cell2mat(sdesc'));
  if iprev
      sdesc_new1(end+1) = sdesc(iprev);
      sdesc_new2(end+1) = sdesc_np(i);
  elseif inext
      sdesc_new1(end+1) = sdesc(inext);
      sdesc_new2(end+1) = sdesc_np(i);
  end
end


sdesc_np = setdiff(sdesc,sdesc_new1);
for i=1:numel(sdesc_np)
      sdesc_new1(end+1) = sdesc_np(i);
      sdesc_new2(end+1) = {repmat('Z',1,lsdesc)};
end
sdesc_np = setdiff(sdesc2,sdesc_new2);
for i=1:numel(sdesc_np)
      sdesc_new1(end+1) = {repmat('Z',1,lsdesc)};
      sdesc_new2(end+1) = sdesc_np(i);
end

N = numel(sdesc_new1);

T(1,1)={'conf1'};
T(1,2)={'conf2'};

%flist=[{'tbeta'} {'tgamma'} {'tdelta'} {'tepsilon'} {'Pdeg'} {'tchi'} {'numax'} {'tau01'} {'tau02'} {'tau03'} {'maxtorinbase'}];
flist=[{'tbeta'} {'tgamma'} {'tdelta'} {'tepsilon'} {'Pdeg'}  {'tchi'} {'numax'} ];
for i=1:numel(flist)
 T(1,3+3*(i-1))={[flist{i} '1']};
 T(1,4+3*(i-1))={[flist{i} '2']};
 T(1,5+3*(i-1))={['abs d' flist{i}]};
end
T(1,3+3*(i))={['dGOenergy' '1']};
T(1,4+3*(i))={['dGOenergy' '2']};
T(1,5+3*(i))={['abs d' 'dGOenergy']};
T(1,3+3*(i+1))={['DM' '1']};
T(1,4+3*(i+1))={['DM' '2']};
T(1,5+3*(i+1))={['abs d' 'DM']};

for i=1:N

  j1=strcmpcellar(sdesc,sdesc_new1{i});
  j2=strcmpcellar(sdesc2,sdesc_new2{i});

  if ~isempty(j1)
      T(i+1,1)={workdb1(j1).prop.sdesc};
  else
      T(i+1,1)={''};
  end
  if ~isempty(j2)
      T(i+1,2)={workdb2(j2).prop.sdesc};
  else
      T(i+1,2)={''};
  end

  for ii=1:numel(flist)
    if ~isempty(j1)
        val1=getfield(workdb1(j1).prop,flist{ii});
    else
        val1=0;
    end
    if ~isempty(j2)
        val2=getfield(workdb2(j2).prop,flist{ii});
    else
        val2=0;
    end
    dval=val1-val2;
    dval=dval-360*fix(dval/180);

    T(i+1,3+3*(ii-1))={val1};
    T(i+1,4+3*(ii-1))={val2};
    T(i+1,5+3*(ii-1))={abs(dval)};

  end

    val1=dGOenergy1(j1);
    val2=dGOenergy2(j2);
    dval=val1-val2;
    T(i+1,3+3*(ii))={val1};
    T(i+1,4+3*(ii))={val2};
    T(i+1,5+3*(ii))={abs(dval)};

    val1=DM1(j1);
    val2=DM2(j2);
    dval=val1-val2;
    T(i+1,3+3*(ii+1))={val1};
    T(i+1,4+3*(ii+1))={val2};
    T(i+1,5+3*(ii+1))={abs(dval)};

end


  
xlswrite(xlsfile,T,'geom comparison','A1');
disp('Done!')
