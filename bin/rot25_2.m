%calculate angles between glycosidic bond and DM
%for r8
%
% Version 1.0    
% Last modified  R O Zhurakivsky 2008-05-19
% Created        R O Zhurakivsky 2006-?-?

atomsind


workdbname=[CD.dbdir filesep 'r8_g.mat']

load(workdbname,'workdb')

recnum=numel(workdb);

[glycvect,DMvect,N1O2vect]=deal(zeros(0));
for i=1:recnum
  if workdb(i).new=='Y'

    ms0=workdb(i);

    iC1 = ms0.ind(find(find(strcmp(pind.labels,'pC1'))==ms0.pind)); 
    iH11= ms0.ind(find(find(strcmp(pind.labels,'pH11'))==ms0.pind)); 

    j=size(glycvect,1)+1;
    sdesc(j,:)=ms0.prop.sdesc;
    glycvect(j,:)=createvect(ms0,iC1,ms0,iH11);
	cos_glycvect(j,:)=glycvect(j,:)/norm(glycvect(j,:));
    DMvect(j,:)=ms0.gaussian.DMv;

if 1
    cla;
	plotmol(ms0);
	plot3(0,0,0,'rs');
	plot3([-DMvect(end,1)/2,DMvect(end,1)/2],[-DMvect(end,2)/2,DMvect(end,2)/2],[-DMvect(end,3)/2,DMvect(end,3)/2],'r');
	plot3(DMvect(end,1)/2,DMvect(end,2)/2,DMvect(end,3)/2,'r*');
	rotate3d on
	angDMglyc=acosd(dot(DMvect(j,:),glycvect(j,:),2)/norm(DMvect(j,:))/norm(glycvect(j,:)));
	[angDMglyc norm(DMvect(j,:))*cosd(angDMglyc)]		
	pause
end

  end
end

%projection of DM vector on C1H11 direction
projDMglyc = dot(DMvect,cos_glycvect,2)


if 0

dglycvect=sqrt(sum(glycvect.^2,2));
dDMvect=sqrt(sum(DMvect.^2,2));

ang1=acosd(dot(glycvect,DMvect,2)./dglycvect./dDMvect)


j1 = glycvect./repmat(dglycvect,1,3);
j3 = cross(N1O2vect,j1);
j3 = j3./repmat(sqrt(sum(j3.^2,2)),1,3);
j2 = cross(j3,j1);


ang2=acosd(dot(DMvect,j2,2)./dDMvect)

bmin=1.1*min([glycvect;DMvect;N1O2vect]);
bmax=1.1*max([glycvect;DMvect;N1O2vect]);
axis([bmin(1) bmax(1) bmin(2) bmax(2) bmin(3) bmax(3)]);

for i=1:size(DMvect,1)
cla
    hold on
    plot3([0 glycvect(i,1)],[0 glycvect(i,2)],[0 glycvect(i,3)],'g');
    text(glycvect(i,1),glycvect(i,2),glycvect(i,3),'C1N1');
    plot3([0 N1O2vect(i,1)],[0 N1O2vect(i,2)],[0 N1O2vect(i,3)],'b');
    text(N1O2vect(i,1),N1O2vect(i,2),N1O2vect(i,3),'N1O2');
    plot3([0 DMvect(i,1)],[0 DMvect(i,2)],[0 DMvect(i,3)],'r');
    text(DMvect(i,1),DMvect(i,2),DMvect(i,3),'DM');
	title([int2str(i) ' : ' workdb(i).prop.sdesc]);
rotate3d on

pause
end

end