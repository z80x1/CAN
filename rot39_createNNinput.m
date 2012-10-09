function res=rot39_createNNinput()
%create input for neuron nets simulations
%
% Version 1.0    
% Last modified  R O Zhurakivsky 2011-10-27
% Created        R O Zhurakivsky 2011-10-27


tic
%clear 
format compact

pindsdef
atomsind

global pind;

%moltype 	= 8 %#ok
%theory 		= 'dftV2'  %#ok
%energy_field = 'MP2_6311__G2dfpd';
%----------------------------------------
moltype 	= 7 %#ok
theory 		= 'dftV1'  %#ok
energy_field = 'MP2_6311__Gdp';

usedpackage = 'Gaussian'  %#ok

onlyoriginal 	= 1;  
outdir = CD.xlsdir;
T = 298.15;
%----------------------------------------


workdbname=['r' int2str(moltype)];
if strcmp(usedpackage,'Gaussian')
  workdbname=[workdbname '_g'];
end
if ~strcmp(theory,'dft')
  workdbname=[workdbname '_' theory];
end
if onlyoriginal
    templ='_or';
    workdbname = [workdbname templ];
end
outfile=[outdir filesep workdbname '_nn.csv']  %#ok
workdbname = [CD.dbdir filesep workdbname '.mat'] %#ok

load(workdbname,'workdb');

recnum=numel(workdb);
istart=1;

clear res;

for iconf=istart:recnum

    mol = workdb(iconf);

    disp([int2str(iconf) ', ' mol.prop.sdesc])

    orderanchor=[4 6 1]; %pinds pC4, pO4, pC1
    order0=[];
    for I=1:numel(orderanchor)
      order0(end+1)=find(mol.pind==orderanchor(I)); %hard indexes
    end
    order=createbondchain(mol,order0,1); %hard indexes   
    
%    if iconf==istart  %zhr111027 не играет роли
%        orderanchor = mol.pind(order);
%    end
    
    mol=createzmt(mol,order);

%    clear s;
    s = struct();
    s.sdesc = mol.prop.sdesc;

%	zmt={};
%    zmt(end+1) = mol.labels(order(1));
%    zmt(end+1:end+3)={mol.labels{order(2)} 1 mol.R(order(2))};
%    i = order(3);
%    zmt(end+1:end+5)={mol.labels{i} 2 mol.R(i) 1 mol.alfa(i)};
%    for i=order(4:end)
%	  zmt(end+1:end+7)={mol.labels{i} mol.iR(i) mol.R(i) mol.ialfa(i) mol.alfa(i) mol.ibeta(i) mol.beta(i)};
%    end
%	res(iconf).zmt = zmt;
    
    s.(['r' an(mol,order(1)) an(mol,order(2))]) = mol.R(order(2));

    s.(['r' an(mol,order(2)) an(mol,order(3))]) = mol.R(order(3));
    s.(['a' an(mol,order(1)) an(mol,order(2)) an(mol,order(3))]) = mol.alfa(order(3));
    
    for i=4:numel(order)
        s.(['r' an(mol,order(i-1)) an(mol,order(i))]) = mol.R(order(i));
        s.(['a' an(mol,order(i-2)) an(mol,order(i-1)) an(mol,order(i))]) = mol.alfa(order(i));
        s.(['d' an(mol,order(i-3)) an(mol,order(i-2)) an(mol,order(i-1)) an(mol,order(i))]) = mol.beta(order(i));
    end
    
    s.DM = mol.gaussian.DM;

    s.E_GO = mol.gaussian.B3LYP_631Gdp; %a.e.

    s.E_SP = getfield(mol.gaussian,energy_field); %a.e.

    indT = find(mol.gaussian.T==T);
	GEC = mol.gaussian.GEC(indT);
    s.GEC = GEC;

    s.G_SP = s.E_SP+GEC/CC.encoef; %a.e.

    res(iconf) = s;
end

E={};
[E{1:length(res),1}] = deal(res.E_SP);
E=cell2mat(E);
dE=(E-min(E))*CC.encoef; %kcal/mol
for j=1:numel(res)
    res(j).dE_SP = dE(j);
end

G={};
[G{1:length(res),1}] = deal(res.G_SP);
G=cell2mat(G);
dG=(G-min(G))*CC.encoef; %kcal/mol
for j=1:numel(res)
    res(j).dG_SP = dG(j);
end

res = rmfield(res,'GEC');

rcelltitles = fieldnames(res);
rcell = struct2cell(res); % 55x1x56
rcell = squeeze(rcell); % 55x56

rcell2 = cell(numel(res)+1,numel(rcelltitles)); %  57x55
rcell2(1,:)=rcelltitles; %  55x1 
rcell2(2:end,:)=rcell'; % 

%res(1,:) = {'sdesc', 'Z-matrix', 'Z-matrix variables only', 'dipole moment,D', ...
%            'E at B3LYP/6-31G(d,p),a.e.', 'E at MP2/6-311++G(d,p)//B3LYP/6-31G(d,p),a.e.', ...
%            'G at MP2/6-311++G(d,p)//B3LYP/6-31G(d,p),a.e.'};
cell2csv(outfile,rcell2,';',2003);       

toc
disp('Done!')
return;

%-------------------------
function name = an(mol, order)
%form atom name

global pind;
name = pind.labels{mol.pind(order)};
if name(1)=='p'
    name=name(2:end);
end    
return;
