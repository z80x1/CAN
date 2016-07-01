function  msfreq=freqrestruct(ms0,moltype);
%sorting atoms in frequency calculation structure in main molecule structure order
%
% Version 1.0    
% Last modified  R O Zhurakivsky 2007-05-01 
% Created        R O Zhurakivsky 2005-09-?

msfreq = ms0.freq;
msfreq.atomnum=numel(msfreq.labels);
msfreq.ind=1:msfreq.atomnum;
msfreq=createbondtable(msfreq);
msfreq=identmol(msfreq,moltype);
[XXX,I1]=sort(ms0.pind);
[XXX,I2]=sort(msfreq.pind);
[X,III]=sortrows([I1,I2]);
I=X(:,2);

%as we obtain I sort vector we have no need in temporary msfreq structure: reusing this name
msfreq = ms0.freq;
msfreq.atomnum=numel(msfreq.labels);
msfreq.ind=1:msfreq.atomnum;
msfreq.labels=msfreq.labels(I);
msfreq.x=msfreq.x(I);
msfreq.y=msfreq.y(I);
msfreq.z=msfreq.z(I);
msfreq.dx=msfreq.dx(:,I);
msfreq.dy=msfreq.dy(:,I);
msfreq.dz=msfreq.dz(:,I);
msfreq=createbondtable(msfreq);
