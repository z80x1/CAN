function ms0=calcproperties(ms0,moltype)
%calculate molecule proporties
%
% Version 1.0    
% Last modified  R O Zhurakivsky 2009-01-13
% Created        R O Zhurakivsky 2005-09-?

%2006-0323 !!changed definition of chi angle from prop.tchi>0 to abs(prop.tchi)<90
%2006-0608 prop.P field excluded
%2006-0718 changed prop.thi field name to prop.tchi
%2006-0811 changed way of confdesc calculations
%2008-0119 changed determination of tchi torsion
%2009-0113 now most of properties are calculated only if parameters they
%depend on exist

if nargin<2
    error('moltype is not specified!');
%    moltype=7;
end

global GLsdesc GLsdescr
global gl_fl_oldchistyle
pindsdef
atomsind

    iC1  = ms0.ind(find(find(strcmp(pind.labels,'pC1'))==ms0.pind)); 
    iC2  = ms0.ind(find(find(strcmp(pind.labels,'pC2'))==ms0.pind)); 
    iC3  = ms0.ind(find(find(strcmp(pind.labels,'pC3'))==ms0.pind)); 
    iC4  = ms0.ind(find(find(strcmp(pind.labels,'pC4'))==ms0.pind)); 
    iC5  = ms0.ind(find(find(strcmp(pind.labels,'pC5'))==ms0.pind)); 
    iO2  = ms0.ind(find(find(strcmp(pind.labels,'pO2'))==ms0.pind)); 
    iO3  = ms0.ind(find(find(strcmp(pind.labels,'pO3'))==ms0.pind)); 
    
    if moltype>100 && mod(moltype,100)==11 %tionucleosides
        iO4  = ms0.ind(find(find(strcmp(pind.labels,'pS4'))==ms0.pind)); 
    else
        iO4  = ms0.ind(find(find(strcmp(pind.labels,'pO4'))==ms0.pind)); 
    end
    
    iO5  = ms0.ind(find(find(strcmp(pind.labels,'pO5'))==ms0.pind)); 
    iH11 = ms0.ind(find(find(strcmp(pind.labels,'pH11'))==ms0.pind)); 
    iH12 = ms0.ind(find(find(strcmp(pind.labels,'pH12'))==ms0.pind)); 
    iH21 = ms0.ind(find(find(strcmp(pind.labels,'pH21'))==ms0.pind)); 
    iH22 = ms0.ind(find(find(strcmp(pind.labels,'pH22'))==ms0.pind)); 
    iH31 = ms0.ind(find(find(strcmp(pind.labels,'pH31'))==ms0.pind)); 
    iH32 = ms0.ind(find(find(strcmp(pind.labels,'pH32'))==ms0.pind)); 
    iH53 = ms0.ind(find(find(strcmp(pind.labels,'pH53'))==ms0.pind)); 

    ibH1 = ms0.ind(find(find(strcmp(pind.labels,'bH1'))==ms0.pind)); 
    ibN1 = ms0.ind(find(find(strcmp(pind.labels,'bN1'))==ms0.pind)); 
    ibC2 = ms0.ind(find(find(strcmp(pind.labels,'bC2'))==ms0.pind)); 
    ibN3 = ms0.ind(find(find(strcmp(pind.labels,'bN3'))==ms0.pind)); 
    ibN  = ms0.ind(find(find(strcmp(pind.labels,'bN'))==ms0.pind));  %N atom of amine group of Cyt
    ibH41= ms0.ind(find(find(strcmp(pind.labels,'bH41'))==ms0.pind)); 
    ibH42= ms0.ind(find(find(strcmp(pind.labels,'bH42'))==ms0.pind)); 

    ibCN2= ms0.ind(find(find(strcmp(pind.labels,'bCN2'))==ms0.pind));  %N atom of amine group of Guo
    ibH21= ms0.ind(find(find(strcmp(pind.labels,'bH21'))==ms0.pind)); 
    ibH22= ms0.ind(find(find(strcmp(pind.labels,'bH22'))==ms0.pind)); 
    ibCN6= ms0.ind(find(find(strcmp(pind.labels,'bCN6'))==ms0.pind));  %N atom of amine group of Ade
    ibH61= ms0.ind(find(find(strcmp(pind.labels,'bH61'))==ms0.pind)); 
    ibH62= ms0.ind(find(find(strcmp(pind.labels,'bH62'))==ms0.pind)); 

    ibC4 = ms0.ind(find(find(strcmp(pind.labels,'bC4'))==ms0.pind)); 
    ibC5 = ms0.ind(find(find(strcmp(pind.labels,'bC5'))==ms0.pind)); 

    ibN7 = ms0.ind(find(find(strcmp(pind.labels,'bN7'))==ms0.pind)); 
    ibC8 = ms0.ind(find(find(strcmp(pind.labels,'bC8'))==ms0.pind)); 
    ibN9 = ms0.ind(find(find(strcmp(pind.labels,'bN9'))==ms0.pind)); 
    ibN8 = ms0.ind(find(find(strcmp(pind.labels,'bN8'))==ms0.pind)); 

    ibF6 = ms0.ind(find(find(strcmp(pind.labels,'bF6'))==ms0.pind)); 
    ibBr5 = ms0.ind(find(find(strcmp(pind.labels,'bBr5'))==ms0.pind)); 

    if any(moltype==[19,21]) %d5AU, d5AC
      ibN5= ms0.ind(find(find(strcmp(pind.labels,'bN5'))==ms0.pind)); 
    else
      ibC5= ms0.ind(find(find(strcmp(pind.labels,'bC5'))==ms0.pind)); 
    end
%    if any(moltype==[11,14]) %only for d6AC & 6AC
      ibN6= ms0.ind(find(find(strcmp(pind.labels,'bN6'))==ms0.pind)); 
%    elseif any(moltype==[9,21,230]) %only for cytidine (dCyt, d5AC)
      ibC6= ms0.ind(find(find(strcmp(pind.labels,'bC6'))==ms0.pind)); 
%    end
    iP5= ms0.ind(find(find(strcmp(pind.labels,'pP5'))==ms0.pind)); 
    iP3= ms0.ind(find(find(strcmp(pind.labels,'pP3'))==ms0.pind)); 
%--------------------------------------------------------------------


      
%    if (moltype>100 && mod(moltype,100)==30) || moltype==22 % 5'nucleotides
    if ~(isempty(iP5) || isempty(iO5) || isempty(iC5) || isempty(iC4))
      prop.tbeta=torang(ms0,iP5,iO5,iC5,iC4);
%    elseif ~any(moltype==[4,5,30])
    elseif ~(isempty(iH53) || isempty(iO5) || isempty(iC5) || isempty(iC4))
      prop.tbeta=torang(ms0,iH53,iO5,iC5,iC4);
    end

%    if ~any(moltype==[4,5,30])
    if ~(isempty(iO5) || isempty(iC5) || isempty(iC4) || isempty(iC3))
      prop.tgamma=torang(ms0,iO5,iC5,iC4,iC3);
    end
     
%    if (all(moltype~=[6,30])) && mod(moltype,100)~=12 %non 5'nucleotides, 3'deoxy or dideoxy molecules 
    if ~(isempty(iO4) || isempty(iC4) || isempty(iC3) || isempty(iO3))
      prop.tdelta_old=torang(ms0,iO4,iC4,iC3,iO3); %definition as in Seanger (Mir 1984, page 29)
    end
    if ~(isempty(iC5) || isempty(iC4) || isempty(iC3) || isempty(iO3))
      prop.tdelta=torang(ms0,iC5,iC4,iC3,iO3);
    end

    if ~(moltype>100 && mod(moltype,100)==12) %non dideoxy molecules
        if (moltype>100 && mod(moltype,100)==31) || moltype==23 %3'nucleotides
          prop.tepsilon=torang(ms0,iC4,iC3,iO3,iP3);
          prop.tkappa2=torang(ms0,iH31,iC3,iO3,iP3);
        elseif all(moltype~=[6,30]) %non 3'deoxy or dideoxy molecules
          prop.tepsilon=torang(ms0,iC4,iC3,iO3,iH32);
          prop.tkappa2=torang(ms0,iH31,iC3,iO3,iH32);
        end
    end
        
    if any(moltype==[7,14]) || mod(moltype,100)==40 %for ribo- molecules
        prop.teta=torang(ms0,iC3,iC2,iO2,iH22);
        prop.tteta=torang(ms0,iO3,iC3,iC2,iO2);
        prop.tkappa1=torang(ms0,iH21,iC2,iO2,iH22);
    end
    
    if any(moltype==[11,9,12,13,14,17,18,19,21]) ||...
            (moltype>=200 && moltype<=299) || (moltype>=400 && moltype<=499) || (moltype>=500 && moltype<=599) %pirimidines
        prop.tchi=torang(ms0,iO4,iC1,ibN1,ibC2);
        prop.rC1N1=adist(ms0,iC1,ibN1);
    elseif any(moltype==[15,16]) || (moltype>=100 && moltype<=199) || (moltype>=300 && moltype<=399) %purines
        prop.tchi=torang(ms0,iO4,iC1,ibN9,ibC4);
        prop.rC1N9=adist(ms0,iC1,ibN9);
    end

    if any(moltype==[11,9,14,30,240]) %NH2 group for cytidine only (not for all Cyd structures!)
        N1=createplane(ms0,ibN,ms0,ibH42,ms0,ibH41);
        N2=createplane(ms0,ibC2,ms0,ibC4,ms0,ibC5);
        prop.tau10=p2pangle(N1,N2)*180/pi*sign(dot(N1,N2)); %angle between plane of NH2 and nucleobase C2C4C5 plane 

        taubuf=abs(torang(ms0,ibC5,ibC4,ibN,ibH41));
        prop.tau11=min(taubuf,180-taubuf); %out of aminobase plane angle of H41 atom 
        taubuf=abs(torang(ms0,ibC5,ibC4,ibN,ibH42));
        prop.tau12=min(taubuf,180-taubuf); %out of aminobase plane angle of H42 atom 
        
    %    cla
    %   plotmol(ms0)
    %   N1=N1/norm(N1)
    %   N2=N2/norm(N2)
    %   T=cross(N1,N2)
    %   plot3([ms0.x(ibN) ms0.x(ibN)+N1(1)],[ms0.y(ibN) ms0.y(ibN)+N1(2)],[ms0.z(ibN) ms0.z(ibN)+N1(3)],'r');
    %   plot3([ms0.x(ibC4) ms0.x(ibC4)+N2(1)],[ms0.y(ibC4) ms0.y(ibC4)+N2(2)],[ms0.z(ibC4) ms0.z(ibC4)+N2(3)],'m');
    %   plot3([0 T(1)],[0 T(2)],[0 T(3)],'m');
    %   plot3(0,0,0,'ks')
    %   rotate3d on
    %   pause
    end
    if moltype==21 %NH2 group for 5AC
        N1=createplane(ms0,ibN,ms0,ibH42,ms0,ibH41);
        N2=createplane(ms0,ibC2,ms0,ibC4,ms0,ibN5);
        prop.tau10=p2pangle(N1,N2)*180/pi*sign(dot(N1,N2)); %angle between plane of NH2 and nucleicbase C2C4C5 plane 

        taubuf=abs(torang(ms0,ibN5,ibC4,ibN,ibH41));
        prop.tau11=min(taubuf,180-taubuf); %out of aminobase plane angle of H41 atom 
        taubuf=abs(torang(ms0,ibN5,ibC4,ibN,ibH42));
        prop.tau12=min(taubuf,180-taubuf); %out of aminobase plane angle of H42 atom 
    end

    if moltype~=30
      prop.tnu0=torang(ms0,iC4,iO4,iC1,iC2);
      prop.tnu1=torang(ms0,iO4,iC1,iC2,iC3);
      prop.tnu2=torang(ms0,iC1,iC2,iC3,iC4);
      prop.tnu3=torang(ms0,iC2,iC3,iC4,iO4);
      prop.tnu4=torang(ms0,iC3,iC4,iO4,iC1);
      prop.rnu0=adist(ms0,iO4,iC1);
      prop.rnu1=adist(ms0,iC1,iC2);
      prop.rnu2=adist(ms0,iC2,iC3);
      prop.rnu3=adist(ms0,iC3,iC4);
      prop.rnu4=adist(ms0,iC4,iO4);
      prop.aO4C1C2=valang(ms0,iO4,iC1,iC2);
      prop.aC1C2C3=valang(ms0,iC1,iC2,iC3);
      prop.aC2C3C4=valang(ms0,iC2,iC3,iC4);
      prop.aC3C4O4=valang(ms0,iC3,iC4,iO4);
      prop.aC4O4C1=valang(ms0,iC4,iO4,iC1);
    
      prop.sumofnu = sum([prop.tnu0,prop.tnu1,prop.tnu2,prop.tnu3,prop.tnu4]);

%2008-0208: corrected Pdeg - 180deg are added when tnu2<0
      %[Saenger,1984 p.31; Altona,Sandaralingam JACS 94,p.8205,1972]
      nom = ((prop.tnu4+prop.tnu1)-(prop.tnu3+prop.tnu0));
      denom = (2*prop.tnu2*(sind(36)+sind(72)));
      P = atan( nom/denom ); 
      if prop.tnu2<0
        P=P+pi;
      end
      prop.Pdeg = mod(P/pi*180,360);
      prop.numax = prop.tnu2/cosd(prop.Pdeg);

%      prop.Pdeg = mod(P/pi*180,180);
%      if prop.numax<0
%      prop.Pdeg=prop.Pdeg+180;
%           prop.numax=-prop.numax;
%      end     

    end

%if moltype==7 | moltype==8
%    prop.O5Hn=adist(ms0,iO5,iH12); %?
%elseif moltype==11 | moltype==9 | moltype==12 | moltype==13 | moltype==14 | moltype==21
%    prop.O5Hn=adist(ms0,iO5,ibN1); %?
%end


    Hbonddata=[];
    Hbonddatadesc={};
    HbondHind=[];
    HbondO1ind=[];
    HbondO2ind=[];


    %AH...B H-bonds searching by geometry criteria 
    %maybe searching at first for H 
    %HbondHind, HbondO1ind, HbondO2ind - hardindexes

    for i=setdiff(1:ms0.atomnum,strcmpcellar(ms0.labels,'H'))

%        ibond2=find(ms0.btB==i);
%        ibond1=find(ms0.btA==i);
        iAlink = [ms0.btA(ms0.btB==i) ms0.btB(ms0.btA==i)];

        for iH = iAlink  % find hydrogens connected to this atom
            if ms0.labels{iH}=='H'
    %          for j=1:ms0.atomnum % try to find third atom for H bond
    %            if j~=i & ms0.labels(j)~='H'
              for j=setdiff(1:ms0.atomnum,strcmpcellar(ms0.labels,'H')) 
                if j~=i
                  HbondABdist=adist(ms0,i,j);
                  HbondHBdist=adist(ms0,iH,j);
                  HbondAHBang =valang(ms0,i,iH,j);
                  if HbondABdist<CG.lcritABdist && HbondHBdist<CG.lcritHBdist && HbondAHBang>CG.lcritAHBang %verify H bond geometry criteria
                    Hbonddatadesc(end+1)={['r' pind.labels{ms0.pind(i)} pind.labels{ms0.pind(j)}]};
                    Hbonddatadesc(end+1)={['a' pind.labels{ms0.pind(i)} pind.labels{ms0.pind(iH)} pind.labels{ms0.pind(j)}]};
                    Hbonddata(end+1)=HbondABdist;
                    Hbonddata(end+1)=HbondAHBang;
                    HbondHind(end+1)=iH;
                    HbondO1ind(end+1)=i;
                    HbondO2ind(end+1)=j;
                  end
                end
              end
            end
        end
    end
    [prop.Hbonddatadesc,ind]=sort(Hbonddatadesc);
    prop.Hbonddata=Hbonddata(ind);

    ind2=ind(1:numel(HbondHind))/2;
    prop.HbondHind=HbondHind(ind2);
    prop.HbondO1ind=HbondO1ind(ind2);
    prop.HbondO2ind=HbondO2ind(ind2);

    %%prop.HbondHind=unique(HbondHind); %unique is for elimination of repeating H atoms involved in several H bonds


    %extracting frequency values of valency and librational modes of H atoms
    if isfield(ms0,'freq') && isfield(ms0.freq,'dx')

      msfreq=freqrestruct(ms0,moltype);
        
      [HbondLMind,HbondVMind]=deal([]);
%      dev=sqrt(msfreq.dx.^2+msfreq.dy.^2+msfreq.dz.^2); %deviations
%      maxdev=max(dev,[],2);

      %searching for librational and stretching(valence) modes of H atoms
      for i=strcmpcellar(msfreq.labels,'H')

%        ibondB=find(msfreq.btB==i);
%        ibondA=find(msfreq.btA==i);
        iOlink = [msfreq.btA(msfreq.btB==i) msfreq.btB(msfreq.btA==i)]; %this may be C atom too

        bonddir=[msfreq.x(i)-msfreq.x(iOlink) msfreq.y(i)-msfreq.y(iOlink) msfreq.z(i)-msfreq.z(iOlink)];
        dr=[msfreq.dx(:,i) msfreq.dy(:,i) msfreq.dz(:,i)];
        bonddir=repmat(bonddir,size(dr,1),1);
        [XXX,I]=max(abs(dot(bonddir,dr,2))); %stretching(valence) mode is determined as mode with maximal deviation in OH (CH) direction
        HbondVMind(end+1)=I;

        S=cross(bonddir,dr,2); %librational mode is determined as mode with maximal deviation perpendicular to OH (CH) direction
        [XXX,I]=max(S(:,1).^2+S(:,2).^2+S(:,3).^2);
        HbondLMind(end+1)=I;
        
    %    [xxx,sortlist]=sort(dev(:,i),1,'descend');
    %    modesind=(sortlist(1:2)); % select two modes with greatest H atom amplitude
    %    HbondLMind(end+1)=min(modesind); %index of libration mode for this H atom
    %    HbondVMind(end+1)=max(modesind); %index of stretching(valence) mode for this H atom
      end
      prop.HbondLMind=HbondLMind;
      prop.HbondVMind=HbondVMind;
%    else
%      prop.HbondLMind=0;
%      prop.HbondVMind=0;
    end

    if any(moltype==[7,14]) %non 2'deoxy molecules only
        prop.rO2H2 = adist(ms0,iO2,iH22);
        prop.rO2O3 = adist(ms0,iO2,iO3);
        prop.aO2H2O3 = valang(ms0,iO2,iH22,iO3);
        prop.rO2O4 = adist(ms0,iO2,iO4);
        prop.aO2H2O4 = valang(ms0,iO2,iH22,iO4);
        prop.rO2O5 = adist(ms0,iO2,iO5);
        prop.aO2H2O5 = valang(ms0,iO2,iH22,iO5);
    end    

    if ~( (moltype>100 && any(mod(moltype,100)==[31,12])) || any(moltype==[6,23,30]) ) %not for 3'nucleotides
        prop.rO3H3 = adist(ms0,iO3,iH32);
        prop.rO3O4 = adist(ms0,iO3,iO4);
        prop.aO3H3O4 = valang(ms0,iO3,iH32,iO4);
    end
%    if ~( (moltype>100 && any(mod(moltype,100)==[31,12])) || any(moltype==[4,5,6,23,30]) ) %not for 3'nucleotides
    if ~(isempty(iO3) || isempty(iO5))
        prop.rO3O5 = adist(ms0,iO3,iO5);
    end
    if ~(isempty(iO3) || isempty(iH32) || isempty(iO5))
        prop.aO3H3O5 = valang(ms0,iO3,iH32,iO5);
    end

    if any(moltype==[7,14]) %non 2'deoxy molecules only
        prop.rO3O2 = adist(ms0,iO3,iO2);
        prop.aO3H3O2 = valang(ms0,iO3,iH32,iO2);
    end

%    if ~((moltype>100 && mod(moltype,100)==30) || any(moltype==[4,5,22,30])) %not for 5'nucleotides
    if ~(isempty(iO5) || isempty(iH53))
      prop.rO5H5 = adist(ms0,iO5,iH53);
    end
%      if ~( moltype==6 || mod(moltype,100)==12 ) %non dideoxy or 3'deoxy molecules
    if ~(isempty(iO5) || isempty(iO3))
          prop.rO5O3 = adist(ms0,iO5,iO3);
    end
    if ~(isempty(iO5) || isempty(iH53) || isempty(iO3))
          prop.aO5H5O3 = valang(ms0,iO5,iH53,iO3);
    end
      
    if ~(isempty(iO5) || isempty(iO4))
      prop.rO5O4 = adist(ms0,iO5,iO4);
    end
    if ~(isempty(iO5) || isempty(iH53) || isempty(iO4))
      prop.aO5H5O4 = valang(ms0,iO5,iH53,iO4);
    end

    if ~any(moltype==30) %not for bases alone
      if ~(isempty(iC2) || isempty(iH21))
          prop.rC2H21 = adist(ms0,iC2,iH21);
      end

      if ~(isempty(iO5) || isempty(iC2))
          prop.rO5C2 = adist(ms0,iO5,iC2);
      end
      if ~(isempty(iO5) || isempty(iC2) || isempty(iH21))
          prop.aO5H21C2 = valang(ms0,iO5,iH21,iC2);
      end
      if ~(isempty(iC3) || isempty(iH31))
          prop.rC3H31 = adist(ms0,iC3,iH31);
      end
      if ~(isempty(iO5) || isempty(iC3))
          prop.rO5C3 = adist(ms0,iO5,iC3);
      end
      if ~(isempty(iO5) || isempty(iC3) || isempty(iH31))
          prop.aO5H31C3 = valang(ms0,iO5,iH31,iC3);
      end
      if ~(isempty(iC1) || isempty(iH11))
        prop.rC2H11 = adist(ms0,iC1,iH11);
      end
      if ~(isempty(iO5) || isempty(iC1))
        prop.rO5C1 = adist(ms0,iO5,iC1);
      end
      if ~(isempty(iO5) || isempty(iC1) || isempty(iH11))
        prop.aO5H11C1 = valang(ms0,iO5,iH11,iC1);
      end
      if ~(isempty(iC2) || isempty(iC1) || isempty(iO4) || isempty(iC4))
          prop.dplC2= dist2plane(ms0,iC2,iC1,iO4,iC4); %distance from plane C1'O4'C4' to atom C2'
      end
      if ~(isempty(iC3) || isempty(iC1) || isempty(iO4) || isempty(iC4))
          prop.dplC3= dist2plane(ms0,iC3,iC1,iO4,iC4); %distance from plane C1'O4'C4' to atom C3'
      end

%build plane over all 5 ring atoms - incorrect
%       aindx=[iC1,iC2,iC3,iC4,iO4]; %indexes of atoms in aminobase plane
%        [x0,a,d,normd] = lsplane([ms0.x(aindx) ms0.y(aindx) ms0.z(aindx)]); %build least square plane throw ALL ring atoms
%       prop.dplC2=d(2); %distance from least square fit ring plane to atom C2'
%       prop.dplC3=d(3); %distance from least square fit ring plane to atom C3'


      if 0
    % determining sugar ring conformation by minimal distances form two atoms to the plane with other three ones.
        ringat = [iC1 iC2 iC3 iC4 iO4];
        for ii=1:5

            ringd(1,ii)=dist2plane(ms0,ringat(4),ringat(1),ringat(2),ringat(3));
            ringd(2,ii)=dist2plane(ms0,ringat(5),ringat(1),ringat(2),ringat(3));
        
            ringat=circshift(ringat,[0 -1]);
        end
        [val,jj]=max(abs(ringd));
        [val,ii]=min(val);
        jj=jj(ii);

        prop.confdesc=pind.labels(find( pind.ind==ms0.pind(ringat(mod(ii+jj-3-1,5)+1)) ));
        prop.confdesc=prop.confdesc{1}(2:end);

        if ringd(jj,ii)>0
          prop.confdesc=strcat(prop.confdesc,'exo');
        else
          prop.confdesc=strcat(prop.confdesc,'endo');
        end
      else
    %determining sugar ring conformation only by P angle value

        confdescs = [{'C3endo'},{'C4exo'},{'O4endo'},{'C1exo'},{'C2endo'},{'C3exo'},{'C4endo'},{'O4exo'},{'C1endo'},{'C2exo'}];

        prop.confdesc = confdescs(ceil(prop.Pdeg/36));
      end


      %???why not using torang function ??
      if ~(isempty(iC5) || isempty(iO5) || isempty(iC4) || isempty(iC3))
          prop.hta2=vecangd(ms0,iC5,iO5,iC4,iC3);  %Hydroden torsion angle  
      end
      if ~(isempty(iO5) || isempty(iP5) || isempty(iC5) || isempty(iC4))
        prop.hta3=vecangd(ms0,iO5,iP5,iC5,iC4);
      elseif ~(isempty(iO5) || isempty(iH53) || isempty(iC5) || isempty(iC4))
        prop.hta3=vecangd(ms0,iO5,iH53,iC5,iC4);
      end
      
      if ~(isempty(iO3) || isempty(iP3) || isempty(iC3) || isempty(iH31))
            prop.hta4=vecangd(ms0,iO3,iP3,iC3,iH31);
      elseif ~(isempty(iO3) || isempty(iH32) || isempty(iC3) || isempty(iH31))
            prop.hta4=vecangd(ms0,iO3,iH32,iC3,iH31);
      end
      
      if ~(isempty(iO2) || isempty(iH22) || isempty(iC2) || isempty(iH21))
        prop.hta5=vecangd(ms0,iO2,iH22,iC2,iH21);
      end
      if ~(isempty(iC3) || isempty(iC2) || isempty(iC4) || isempty(iC5))
        prop.hta0=vecangd(ms0,iC3,iC2,iC4,iC5);  % characterizes difference between alfa & beta form of ribose
      end
    end

%    if any(moltype==[11,9,14,21,12,13,17,18,19,30,  15,16]) || any(floor(moltype/100)==[2,4,5])
    %only for cytidine, thymidine, uridine
    if ~any(moltype==[2,3,4,5,6,7,8]) %now this block is for all

        if moltype<100
            if any(moltype==[9,12,13]) %dCyd, dUrd, dThd
                aindx=[ibN1,ibC2,ibN3,ibC4,ibC5,ibC6]; %indexes of atoms in aminobase plane
            elseif any(moltype==[19,21]) %d5AU , d5AC
                aindx=[ibN1,ibC2,ibN3,ibC4,ibN5,ibC6];
            elseif any(moltype==[11,14,17,18,30]) % 6Ax
                aindx=[ibN1,ibC2,ibN3,ibC4,ibC5,ibN6];
            elseif  any(moltype==[15,16]) % dAde, dGuo
                aindx =[ibN1,ibC2,ibN3,ibC4,ibC5,ibC6];
                aindx2=[ibC4,ibC5,ibN7,ibC8,ibN9];
            else
                error('errorD001: Type for aindx is not specified correctly');
            end
        else %moltype >=100
            if any(floor(moltype/100)==[2,4,5]) %Cyd, Urd, Thd
                if any(mod(moltype,100)==[30,31,40,50,51,52,61,70,11,12,21,22]) 
                    aindx=[ibN1,ibC2,ibN3,ibC4,ibC5,ibC6]; %indexes of atoms in aminobase plane
                elseif mod(moltype,100)==20 % 6Fx
                    aindx=[ibN1,ibC2,ibN3,ibC4,ibC5,ibF6];
%                elseif any(mod(moltype,100)==[21,22]) % 5Brx - Br atom is not in ring!!!
%                    aindx=[ibN1,ibC2,ibN3,ibC4,ibBr5,ibC6];
                else
                    error('errorD003: Type for aindx is not specified correctly');
                end
            elseif any(floor(moltype/100)==[1,3]) %Ade, Guo
                if any(mod(moltype,100)==[30,31,40,50,51,52,53,61,70,11,12,90,91])
                    aindx =[ibN1,ibC2,ibN3,ibC4,ibC5,ibC6];
                    aindx2=[ibC4,ibC5,ibN7,ibC8,ibN9];
                elseif any(mod(moltype,100)==[10,62]) %8-aza
                    aindx =[ibN1,ibC2,ibN3,ibC4,ibC5,ibC6];
                    aindx2=[ibC4,ibC5,ibN7,ibN8,ibN9];
                else
                    error('errorD004: Type for aindx is not specified correctly');
                end
            else
                error('errorD002: Type for aindx is not specified correctly');
            end
        end

        
    %    N1bplane=createplane(ms0,ibN1,ms0,ibN3,ms0,ibC5);
    %N1bplane - least square plane of purine ring
        [x0,N1bplane,d,normd]=lsplane([ms0.x(aindx) ms0.y(aindx) ms0.z(aindx)]);

        if any(moltype==30)
          N2=createvect(ms0,ibH1,ms0,ibN1);
        elseif any(moltype==[11,9,12,13,14,17,18,19,21]) ||...
            (moltype>=200 && moltype<=299) || (moltype>=400 && moltype<=499) || (moltype>=500 && moltype<=599) %pirimidines
          N2=createvect(ms0,iC1,ms0,ibN1);
        elseif any(moltype==[15,16]) || (moltype>=100 && moltype<=199) || (moltype>=300 && moltype<=399) %purines
          N2=createvect(ms0,iC1,ms0,ibN9);
        else
          error('Cannot locate atoms of glycosidic bond');
        end

        %%?? maybe this need to be updated to N1 vector obtained by lsplane 
        %sgn=sign(dot(cross(createvect(ms0,aindx(1),ms0,aindx(3)),createvect(ms0,aindx(1),ms0,aindx(5))),N2)); %distinguish cases when bond is under plane or above it
        sgn=sign(dot(N1bplane,N2)); %distinguish cases when bond is under plane or above it
        A=p2pangle(N1bplane,N2)*180/pi*sgn+90;

        %different variants of C1' out of base plane angle
        prop.tau01=A-180*round(A/180); %angle between C1'N1 (or C1'N9) and aminobase least square fit plane  (aglycoutplane)

        %2007-0810 %distance from least square plane of nucleobase to atom C1'
        v1=cross(N1bplane,[1 0 0 ]);
        v2=cross(N1bplane,v1);
        x0=x0';
        x1=x0+v1;
        x2=x0+v2;
        prop.dplC1base= dist2planecart([ms0.x(iC1) ms0.y(iC1) ms0.z(iC1)],x0,x1,x2); 
        
        
        %not only for pirimidines
        prop.tau04=valang(ms0,iC1,aindx(1),aindx(2))+valang(ms0,aindx(2),aindx(1),aindx(6))+valang(ms0,aindx(6),aindx(1),iC1);  %sum of valence angles at N1 atom
        prop.tau05=torang(ms0,iC1,aindx(1),aindx(2),aindx(3));  %torsion angle between sugar and one base side C1'N1C2N3
        prop.tau06=torang(ms0,iC1,aindx(1),aindx(6),aindx(5));  %torsion angle between sugar and another base side C1'N1C6C5


        if any(moltype==[11,9,14,21,30,240]) % for Cyt only

      %    N1bplane=createplane(ms0,ibN1,ms0,ibN3,ms0,ibC5);
      %    [x0,N1bplane,d,normd]=lsplane([ms0.x(aindx) ms0.y(aindx) ms0.z(aindx)]);
          N2=createvect(ms0,ibC4,ms0,ibN);
          %?? maybe this need to be updated to N1 vector obtained by lsplane 
          sgn=sign(dot(cross(createvect(ms0,aindx(1),ms0,aindx(3)),createvect(ms0,aindx(1),ms0,aindx(5))),N2)); %distinguish cases when bond is under plane or above it
          A=p2pangle(N1bplane,N2)*180/pi*sgn+90;
          prop.tau02=A-180*round(A/180); %angle between C4N (or C2N, or C6N) and aminobase least square fit plane

          N1=createplane(ms0,ibN,ms0,ibH41,ms0,ibH42);
          N2=createvect(ms0,ibC4,ms0,ibN);
          sgn=sign(dot(cross(createvect(ms0,ibN,ms0,ibH41),createvect(ms0,ibN,ms0,ibH42)),N2)); %distinguish cases when bond is under plane or above it
          A=p2pangle(N1,N2)*180/pi*sgn+90;
          prop.tau03=A-180*round(A/180); %angle between C4N (or C2N, or C6N) and aminogroup plane NH2
          prop.rCNamino=adist(ms0,ibC4,ibN);
 
        elseif any(moltype==[16,340]) % for Guo only
        
          N2=createvect(ms0,ibC2,ms0,ibCN2);
%          sgn=sign(dot(cross(createvect(ms0,aindx(1),ms0,aindx(3)),createvect(ms0,aindx(1),ms0,aindx(5))),N2)); %distinguish cases when bond is under plane or above it
          sgn=sign(dot(N1bplane,N2)); %distinguish cases when bond is under plane or above it
          A=p2pangle(N1bplane,N2)*180/pi*sgn+90;
          prop.tau02=A-180*round(A/180); %angle between C2N2 and aminobase least square fit plane

          N1=createplane(ms0,ibCN2,ms0,ibH21,ms0,ibH22);
          N2=createvect(ms0,ibC2,ms0,ibCN2);
          sgn=sign(dot(cross(createvect(ms0,ibCN2,ms0,ibH21),createvect(ms0,ibCN2,ms0,ibH22)),N2)); %distinguish cases when bond is under plane or above it
          A=p2pangle(N1,N2)*180/pi*sgn+90;
          prop.tau03=A-180*round(A/180); %angle between C2N2 and aminogroup plane NH2
          prop.rCNamino=adist(ms0,ibC2,ibCN2);

        elseif any(moltype==[15,140]) % for Ado only
          N2=createvect(ms0,ibC6,ms0,ibCN6);
          sgn=sign(dot(N1bplane,N2)); %distinguish cases when bond is under plane or above it
          A=p2pangle(N1bplane,N2)*180/pi*sgn+90;
          prop.tau02=A-180*round(A/180); %angle between C6N6 and aminobase least square fit plane

          N1=createplane(ms0,ibCN6,ms0,ibH61,ms0,ibH62);
          N2=createvect(ms0,ibC6,ms0,ibCN6);
          sgn=sign(dot(cross(createvect(ms0,ibCN6,ms0,ibH61),createvect(ms0,ibCN6,ms0,ibH62)),N2)); %distinguish cases when bond is under plane or above it
          A=p2pangle(N1,N2)*180/pi*sgn+90;
          prop.tau03=A-180*round(A/180); %angle between C6N6 and aminogroup plane NH2
          prop.rCNamino=adist(ms0,ibC6,ibCN6);
        end


        iend = numel(aindx);
        torbpl=zeros(1,iend);
        for i=1:iend
            %for pirimidines for i=1: N6N1C2N3, i=6: C5N6N1C2 
            i1=i-1; if i1==0, i1=iend; end;
            i3=i+1; if i3>iend, i3=i3-iend; end;
            i4=i+2; if i4>iend, i4=i4-iend; end;
            torbpl(i)=torang(ms0,aindx(i1),aindx(i),aindx(i3),aindx(i4)); %torsion angle  in aminobase
        end
%        torbpl(1)=torang(ms0,aindx(6),aindx(1),aindx(2),aindx(3)); %torsion angle N6N1C2N3 in aminobase
        prop.torbpl=torbpl;
        prop.torbplind=aindx;


        if any(moltype==[15,16]) || any(floor(moltype/100)==[1,3]) %Ade, Guo
            iend = numel(aindx2);
            torbpl2=zeros(1,iend);
            for i=1:iend
                i1=i-1; if i1==0, i1=iend; end;
                i3=i+1; if i3>iend, i3=i3-iend; end;
                i4=i+2; if i4>iend, i4=i4-iend; end;
                torbpl2(i)=torang(ms0,aindx2(i1),aindx2(i),aindx2(i3),aindx2(i4)); %torsion angle  in aminobase
            end
            prop.torbpl2=torbpl2;
            prop.torbpl2ind=aindx2;
            prop.maxtorinbase=max(min(abs([torbpl,torbpl2]),180-abs([torbpl,torbpl2])));

        else
            prop.maxtorinbase=max(min(abs(torbpl),180-abs(torbpl))); %maximum torsion angle along the ring(s)

        end

    end

    if ~any(moltype==30)    
      prop.sdesc=GLsdescr(ceil(prop.Pdeg/36));
      if isfield(prop,'hta2')
          prop.sdesc=[prop.sdesc GLsdesc(ceil(prop.hta2/120))]; %gamma:   a (0;120] - g+; b (120;240] - t;  c (240;360] - g-;
      end
      if isfield(prop,'hta3')
          prop.sdesc=[prop.sdesc GLsdesc(ceil(prop.hta3/120))]; %beta:    a (0;120] - g+; b (120;240] - t;  c (240;360] - g-;
      end
      if isfield(prop,'hta4')
          prop.sdesc=[prop.sdesc GLsdesc(ceil(prop.hta4/120))]; %epsilon: a (0;120] - g-; b (120;240] - g+; c (240;360] - t;
      end
      if isfield(prop,'hta5')
        prop.sdesc=[prop.sdesc GLsdesc(ceil(prop.hta5/120))];
      end

      if ~(any(moltype==[2,3,4,5,6,7,8,22,23]))
%       if exist('gl_fl_oldchistyle','var') && gl_fl_oldchistyle
       if ~isempty(gl_fl_oldchistyle) && gl_fl_oldchistyle %gl_fl_oldchistyle is already defined as global

        if abs(prop.tchi)<90
          prop.sdesc=[prop.sdesc 'S'];
        else
          prop.sdesc=[prop.sdesc 'A'];
        end
       else
        %zhr080119
        %as in Saenger (p.33-35)    60deg<chi<=120deg - high-syn (V),  -60deg<chi<=60deg - syn (S)
        %               -150deg<chi<=-60deg - high-anti (B),  -180deg<chi<=-150deg,120deg<chi<=180deg - anti (B)
        if (prop.tchi > -60) && (prop.tchi <= 60)
          prop.sdesc=[prop.sdesc 'T'];
        elseif (prop.tchi > 60) && (prop.tchi <= 120)
          prop.sdesc=[prop.sdesc 'V'];
        elseif (prop.tchi > -150) && (prop.tchi <= -60)
          prop.sdesc=[prop.sdesc 'C'];
        else
          prop.sdesc=[prop.sdesc 'B'];
        end
       end
      end
      
      prop.sdescnew=GLsdescr(ceil(prop.Pdeg/36));
      if isfield(prop,'hta2')
          prop.sdescnew=[prop.sdescnew GLsdesc(ceil(prop.hta2/40))];
      end
      if isfield(prop,'hta3')
          prop.sdescnew=[prop.sdescnew GLsdesc(ceil(prop.hta3/40))];
      end
      if isfield(prop,'hta4')
          prop.sdescnew=[prop.sdescnew GLsdesc(ceil(prop.hta4/40))];
      end
      if isfield(prop,'hta5')
        prop.sdescnew=[prop.sdescnew GLsdesc(ceil(prop.hta5/40))];
      end
    else
      prop.sdesc = 'A';
      prop.sdescnew = 'A';
    end
        
    %create distance and valence angle matrixes for whole molecule
    
    distmatr=zeros(numel(ms0.btA),1);
    distmatrlabels = cell(numel(ms0.btA),1);
    distmatrpinds = zeros(numel(ms0.btA),2,'uint16');
    
    bondorder=1:numel(ms0.btA);
    for i=bondorder
        distmatr(i) = adist(ms0,ms0.ind(ms0.btA(i)),ms0.ind(ms0.btB(i)));
        distmatrpinds(i,:) = sort([ms0.pind(ms0.btA(i)) ms0.pind(ms0.btB(i))]);
%        distmatrlabels(i) = {[pind.labels{ms0.pind(ms0.btA(i))} pind.labels{ms0.pind(ms0.btB(i))}]};
    end

    %creating table that describes all existed valence angles
    %and matrix of torsion angles too
    [angtA,angtB,angtC]=deal(zeros(3*3*numel(ms0.ind),1,'uint16')); %3*3 is creating vector that can contain all possible valence angle of molecule
    [tortA,tortB,tortC,tortD]=deal(zeros(10*numel(ms0.ind),1,'uint16')); %10 is creating vector that can contain all possible torsion angles of molecule
    angtnumel=0;
    tortnumel=0;
    for i=1:numel(ms0.ind)
        angAbonds=[ms0.btB(ms0.btA==i) ms0.btA(ms0.btB==i)];
        for j=angAbonds
            angBbonds=[ms0.btB(ms0.btA==j) ms0.btA(ms0.btB==j)];
            angBbonds=setdiff(angBbonds,i);
            for k=angBbonds
                angtnumel=angtnumel+1;
                angtA(angtnumel)=i;
                angtB(angtnumel)=j;
                angtC(angtnumel)=k;

                torCbonds=[ms0.btB(ms0.btA==k) ms0.btA(ms0.btB==k)];
                torCbonds=setdiff(torCbonds,j);
                for l=torCbonds
                    tortnumel=tortnumel+1;
                    tortA(tortnumel)=i;
                    tortB(tortnumel)=j;
                    tortC(tortnumel)=k;
                    tortD(tortnumel)=l;
                end

            end
        end
    end

    angtA=angtA(1:angtnumel);
    angtB=angtB(1:angtnumel);
    angtC=angtC(1:angtnumel);
    angtAC=sort([angtA angtC],2);
    angt=unique([angtAC(:,1) angtB angtAC(:,2)],'rows');

    anglematr=zeros(size(angt,1),1);
    anglematrlabels = cell(size(angt,1),1);
    anglematrpinds = zeros(size(angt,1),3);
    for i=1:numel(anglematr)
        anglematr(i) = valang(ms0,ms0.ind(angt(i,1)),ms0.ind(angt(i,2)),ms0.ind(angt(i,3)));
        tmp=sort([ms0.pind(angt(i,1)) ms0.pind(angt(i,3))]);
        anglematrpinds(i,:) = [tmp(1) ms0.pind(angt(i,2)) tmp(2)];
%        anglematrlabels(i) = {[pind.labels{ms0.pind(angt(i,1))} pind.labels{ms0.pind(angt(i,2))} pind.labels{ms0.pind(angt(i,3))}]};
    end

    tortA=tortA(1:tortnumel);
    tortB=tortB(1:tortnumel);
    tortC=tortC(1:tortnumel);
    tortD=tortD(1:tortnumel);
    [tortAD,iAD]=sort([tortA tortD],2);
%    tortAD=tortAD';
    tortBC=[tortB tortC];
    for j = 1:tortnumel, tortBC(j,:) = tortBC(j,iAD(j,:)); end
    tort=unique([tortAD(:,1) tortBC tortAD(:,2)],'rows');
    
    tormatr=zeros(size(tort,1),1);
    tormatrlabels = cell(size(tort,1),1);
    tormatrpinds = zeros(size(tort,1),4);
    for i=1:numel(tormatr)
        tormatr(i) = torang(ms0,ms0.ind(tort(i,1)),ms0.ind(tort(i,2)),ms0.ind(tort(i,3)),ms0.ind(tort(i,4)));
        [tmp,itmp]=sort([ms0.pind(tort(i,1)) ms0.pind(tort(i,4))]);
        if itmp(1)==1
            tormatrpinds(i,:) = [tmp(1) ms0.pind(tort(i,2)) ms0.pind(tort(i,3)) tmp(2)];
        else
            tormatrpinds(i,:) = [tmp(1) ms0.pind(tort(i,3)) ms0.pind(tort(i,2)) tmp(2)];
        end
    end

    prop.distmatr = distmatr;
    prop.distmatrpinds = distmatrpinds;
    prop.anglematr = anglematr;
    prop.anglematrpinds = anglematrpinds;
    prop.tormatr = tormatr;
    prop.tormatrpinds = tormatrpinds;

    ms0.prop = prop;
