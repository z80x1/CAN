function ind = findallrotates(ms0,i2,i3,order)
%returns all hardindexes of atoms that are able to rotating around link i2-i3
%order - first indexes of atoms in chain
%
% Version 1.0    
% Last modified  R O Zhurakivsky 2007-11-02 
% Created        R O Zhurakivsky 2005-08-06

%20050806 torI = [ms0.pind... -> torI = [ms0.ind... 
%20070329 added order parameter
%20070719 torI and II all contained now pinds
%20071102 added transformation of ialfa, ibeta, iR with or array

if nargin<4
  error('Incorrect number of parameters used!');
end

or=createbondchain(ms0,order);


%070731 II=sort(ms0.pind(findtorlinks(ms0,i2,i3)),2);
%for r270 new structures created pinds are useless (r27004 creation)
II=sort(ms0.ind(findtorlinks(ms0,i2,i3)),2);
iR=ms0.iR; 		 iR(iR==0)=1;
ialfa=ms0.ialfa; ialfa(ialfa==0)=1;
ibeta=ms0.ibeta; ibeta(ibeta==0)=1;


%061004 PP=sort(ms0.pind); %very bad construction of ZMT matrix
%PP=ms0.pind; %very bad construction of ZMT matrix

%070731 torI = sort([ms0.pind iR ialfa ibeta],2); %for how indexes in ZMT fields are arranged see createbondchain
%for r270 new structures created pinds are useless (r27004 creation)

%torI = sort([uint16(ms0.ind) iR ialfa ibeta],2); %for how indexes in ZMT fields are arranged see createbondchain
torI = sort([uint16(ms0.ind) or(iR)' or(ialfa)' or(ibeta)'],2); %for how indexes in ZMT fields are arranged see createbondchain

[C,IA,ind]=intersect(uint16(II),torI,'rows'); %#ok
