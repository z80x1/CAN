%
% Version 1.0    
% Last modified  R O Zhurakivsky 2006-09-19
% Created        R O Zhurakivsky 2006-?-?



        N1=createplane(ms0,ibN,ms0,ibH42,ms0,ibH41);
        N2=createplane(ms0,ibC2,ms0,ibC4,ms0,ibC5);
        prop.tau10=p2pangle(N1,N2)*180/pi*sign(dot(N1,N2)); %angle between plane of NH2 and nucleicbase C2C4C5 plane 

        taubuf=abs(torang(ms0,ibC5,ibC4,ibN,ibH41));
        prop.tau11=min(taubuf,180-taubuf); %out of aminobase plane angle of H41 atom 
        taubuf=abs(torang(ms0,ibC5,ibC4,ibN,ibH42));
        prop.tau12=min(taubuf,180-taubuf); %out of aminobase plane angle of H42 atom 
      [x0,N1,d,normd]=lsplane([ms0.x(aindx) ms0.y(aindx) ms0.z(aindx)]);
      N2=createvect(ms0,iC1,ms0,ibN1);
      sgn=sign(dot(cross(createvect(ms0,aindx(1),ms0,aindx(3)),createvect(ms0,aindx(1),ms0,aindx(5))),N2)); %distinguish cases when bond is under plane or above it
      A=p2pangle(N1,N2)*180/pi*sgn+90;
      prop.tau01=A-180*round(A/180); %angle between C1'N1 and aminobase least square fit plane  (aglycoutplane)

      if any(moltype==[11,9,14,21]) % for Cyt only
    %    N1=createplane(ms0,ibN1,ms0,ibN3,ms0,ibC5);
    %    [x0,N1,d,normd]=lsplane([ms0.x(aindx) ms0.y(aindx) ms0.z(aindx)]);
        N2=createvect(ms0,ibC4,ms0,ibN);
        sgn=sign(dot(cross(createvect(ms0,aindx(1),ms0,aindx(3)),createvect(ms0,aindx(1),ms0,aindx(5))),N2)); %distinguish cases when bond is under plane or above it
        A=p2pangle(N1,N2)*180/pi*sgn+90;
        prop.tau02=A-180*round(A/180); %angle between C4N and aminobase least square fit plane

        N1=createplane(ms0,ibN,ms0,ibH41,ms0,ibH42);
        N2=createvect(ms0,ibC4,ms0,ibN);
        sgn=sign(dot(cross(createvect(ms0,ibN,ms0,ibH41),createvect(ms0,ibN,ms0,ibH42)),N2)); %distinguish cases when bond is under plane or above it
        A=p2pangle(N1,N2)*180/pi*sgn+90;
        prop.tau03=A-180*round(A/180); %angle between C4N and aminogroup plane NH2
        prop.rCNamino=adist(ms0,ibC4,ibN);
      end
