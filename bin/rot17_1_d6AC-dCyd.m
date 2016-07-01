%convert dAC6 unique conformations to dCyd  ones by changing N6->C6 and adding H6 at 1.084A from C6
%see rot17_5 for new version
%
% Version 1.0    
% Last modified  R O Zhurakivsky 2008-05-19
% Created        R O Zhurakivsky 2006-?-?

clear 
atomsind

workdbname=[CD.dbdir filesep 'r11_g.mat']
gtemplname=[CD.templatesdir filesep 'r901_templ.gjf']
odir=[CD.xyzdir filesep 'xyz'];

global pind
atomsind
load(workdbname,'workdb')

recnum=numel(workdb);
for i=1:recnum

    if workdb(i).new~='Y'
       continue
    end  

    ms0.labels=workdb(i).labels;
    ms0.x=workdb(i).x;
    ms0.y=workdb(i).y;
    ms0.z=workdb(i).z;
    ms0.atomnum=workdb(i).atomnum;
    ms0.ind=workdb(i).ind;
    ms0.desc=workdb(i).prop.sdesc;
    ms0.pind=workdb(i).pind;

    ind=find(ms0.pind==find(strcmp(pind.labels,'bN6')));
    ms0.pind(ind) = find(strcmp(pind.labels,'bC6'));
    ms0.labels(ind)='C';

    indA=find(ms0.pind==find(strcmp(pind.labels,'bN1')));
    indB=ind;
    indC=find(ms0.pind==find(strcmp(pind.labels,'bC5')));
    A=[ms0.x(indA) ms0.y(indA) ms0.z(indA)];
    B=[ms0.x(indB) ms0.y(indB) ms0.z(indB)];
    C=[ms0.x(indC) ms0.y(indC) ms0.z(indC)];
    X=(A+C)./2;

    ortXB=(B-X)/norm(B-X);
    
    dist=1.084;
    D=B+dist*ortXB;

    ms0.atomnum=ms0.atomnum+1;
    ms0.x(end+1)=D(1);
    ms0.y(end+1)=D(2);
    ms0.z(end+1)=D(3);
    ms0.labels(end+1)='H';
    ms0.ind(end+1)=max(ms0.ind)+1;
    ms0.pind(end+1)=find(strcmp(pind.labels,'bH6'));
    ms0=createbondtable(ms0);

    order=[6 1 2]; 

    or = [order, setdiff(1:ms0.atomnum,order)]; %!!if using see rot14 for better inplementation
    savemol(odir,ms0,0,or);

%not needed
%    [xxx,order]=sortrows(ms0.pind);

%   savemolgs(odir,ms0,3,order,gtemplname); %Gaussian with XYZ
    savemolgs(odir,ms0,4,order,gtemplname); %Gaussian with ZMT

end


