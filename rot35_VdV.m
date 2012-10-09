%calculate Van der Vaals volules over conformations for different VdV radii systems
%
% Version 1.0    
% Last modified  R O Zhurakivsky 2008-05-19
% Created        R O Zhurakivsky 2006-?-?

clear

atomsind

format compact

%-------------------------------------
  moltype=9  %#ok
  theory='dftV2'  %#ok
  usedpackage='Gaussian'  %#ok
  onlyoriginal=1;  % process db with only original conformations

%-------------------------------------

workdbname=['r' int2str(moltype)]   %#ok
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
workdbname=[CD.dbdir filesep workdbname '.mat'] %#ok

load(workdbname,'workdb')
recnum=numel(workdb);
if ~recnum
  error('Database is empty')
end



%VDV=VDV.Bondi; disp('Bondi VdV system used');
VDV=VDV.Zorkiy; disp('Zorkiy VdV system used');

VDV_vol=zeros(0);
desc='';

for i=1:recnum
    if workdb(i).new~='Y'
	continue
    end
    ms0=workdb(i);

    desc(end+1,:) = ms0.prop.sdesc;

    V=0;
    for ind=1:numel(ms0.labels)
       V=V+4/3*pi*getfield(VDV,ms0.labels{ind})^3;
    end
    for ind=1:numel(ms0.btA)
	Ri  = getfield(VDV,ms0.labels{ms0.btA(ind)});
	Rj  = getfield(VDV,ms0.labels{ms0.btB(ind)});
	dij = adist(ms0,ms0.btA(ind),ms0.btB(ind));


        hij = Ri - (Ri^2+dij^2-Rj^2)/2/dij;
	dV1=pi*hij^2/3*(3*Ri-hij);

        hji = Rj - (Rj^2+dij^2-Ri^2)/2/dij;
	dV2=pi*hji^2/3*(3*Rj-hji);

	V=V-dV1-dV2;
    end


    VDV_vol(end+1)=V;

%break
end



VDV_volmin=min(VDV_vol) %#ok
VDV_volmax=max(VDV_vol) %#ok
VDV_volmean=mean(VDV_vol) %#ok
VDV_volstd=std(VDV_vol,1) %#ok


disp('Done!')
