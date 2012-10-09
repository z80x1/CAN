%analyze vibrational deviations of atoms
%
% Version 1.0    
% Last modified  R O Zhurakivsky 2008-05-19
% Created        R O Zhurakivsky 2006-?-?

workdbname='r8_g.mat'

load(workdbname,'workdb')

i=1;

ms0=workdb(1);

dev=sqrt(ms0.freq.dx.^2+ms0.freq.dy.^2+ms0.freq.dz.^2); %deviations
maxdev=max(dev,[],2);
dev(dev<0.1*repmat(maxdev,1,ms0.atomnum))=0;   %?

for i=1:numel(ms0.prop.HbondHind)
    Hatomind = ms0.prop.HbondHind(i); % hydrogen atom index in current H bond
    [xxx,sortlist]=sort(dev(:,Hatomind),1,'descend');
    modesind=(sortlist(1:2)) % two modes with greatest H atom amplitude

    ms0.prop.HbondLMind=min(modesind); %index of libration mode for this H atom
    ms0.prop.HbondVMind=max(modesind); %index of stretching(valence) mode for this H atom
end



if 0
dev=sqrt(ms0.freq.dx.^2+ms0.freq.dy.^2+ms0.freq.dz.^2); %deviations
maxdev=max(dev,[],2);
atoms = zeros(size(dev)); %atomsinmotions matrix
atoms(dev>0.5*repmat(maxdev,1,ms0.atomnum))=1;

for i=1:numel(ms0.prop.HbondHind)

    Hatomind = ms0.prop.HbondHind(i); % hydrogen atom index in current H bond
    modes = find(atoms(:,Hatomind));  % modes this H atom is included
    modes1=find(sum(atoms,2)==1); % modes with only one atom vibrating
    ourmodes=intersect(modes,modes1) %modes we are interested
    if numel(ourmodes)~=2
       warning('modes number not equal to 2'); 
    end
    ms0.prop.HbondLMind=ourmodes(1); %index of libration mode for this H atom
    ms0.prop.HbondVMind=ourmodes(2); %index of stretching(valence) mode for this H atom

%    ms0.freq.freq(ourmodes(1))

end
end
