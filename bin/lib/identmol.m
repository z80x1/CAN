function [ms0,status] = identmol(ms0,moltype,fl_restart)
%identifies RIBOSE molecule structure. Returns modified struct ms0 with pind 
%field that contains indexes of atoms as are specified in global structure pind
%also returns table with all molecule bonds ms0.bt
%
% Version 1.0    
% Last modified  R O Zhurakivsky 2009-05-05
% Created        R O Zhurakivsky 2005-09-?

%060517 H11' and H12' are interchanged for non sugar molecules
%060529 H51' and H52' are now distinguished by cross(C5'O5',C5'C4') and C5'H51' angle
%060622 H51' & H52' labels are exchanged to be in agreement with Seanger, H21'-H22' & H31'-H32' label pairs are in disagreement with Seanger
%060806 fl_restart : flag to continue identification
%060812 Thd CH3 hydrogens are distinguished by C4C5CC5Hx torsion values
%060812 code cleaning to MatLAB code checker standarts is performed
%090113 added moltypes 2,3,4,5
%090113 added moltypes 390,391,352,353
%090701 added moltypes 910

%moltype=2: 2'-hydroxy-tetrahydrofurane         (r2)  allien
%moltype=3: 2',3'-hydroxy-tetrahydrofurane      (r3)  allien
%moltype=4: 4'-methyl,2'-hydroxy-tetrahydrofurane    (r4)  allien
%moltype=5: 4'-methyl,2',3'-hydroxy-tetrahydrofurane (r5)  allien
%moltype=6: 1'3'deoxyribose                     (r6)  allien
%moltype=7: 1'deoxyribose                       (r7)  done
%moltype=8: 1'2'deoxyribose                     (r8)  done
%moltype=9: 2'deoxy-cytidine        dCyd        (r9)  done
%moltype=11: 2'deoxy-6azacytidine   d6AC        (rb)  done
%moltype=12: 2'deoxy-uridine        dUrd        (rc)  done
%moltype=13: 2'deoxy-thymidine      dThd        (rd)  done
%moltype=14: 6azacytidine           6AC         (re)
%moltype=15: 2'deoxy-adenosine          dAdo    (rf)  done
%moltype=16: 2'deoxy-guanosine          dGuo    (rg)  done

%moltype=17: 2'deoxy-6azauridine    d6AU        (rh)    2do
%moltype=18: 2'deoxy-6azathymidine  d6AT        (ri)

%moltype=19: 2'deoxy-5azauridine    d5AU        (rj)    2do
%moltype=21: 2'deoxy-5azacytidine   d5AC        (rl)    2do

%moltype=22: 1'2'deoxyribose with P residue at 5' (rm)
%moltype=23: 1'2'deoxyribose with P residue at 3' (rn)

%moltype=24: 1'2'deoxyribose with NH2 residue   (r020)
%moltype=25: 1'2'deoxyribose with OH residue    (r021)

%moltype=30: 6azacytidinebase                   (r030)

%moltype=101: Adenine                           (r101)
%moltype=110: 2'deoxy-8aza-adenosine     d8AA    (r110)
%moltype=111: 2'deoxy-5'tioadenosine            (r111)
%moltype=112: 2',3'-dideoxy-didehydro-adenosine           (r112)
%moltype=130: dAMP 5'                           (r130)
%moltype=131: dAMP 3'                           (r131)
%moltype=140: rAdo                              (r140)
%moltype=150: dAdo imino                        (r150)
%moltype=151: 2'-deoxy-2aminopurine            (r151)       2do
%moltype=152: 2'deoxy-8oxo-adenosine     d8OA   (r152)
%moltype=161: 2'deoxy-purine             dP     (r161)       2do
%moltype=162: 2'-deoxy-8-aza-purine      d8AP   (r162)       2do


%moltype=201: Cytidine                          (r201)
%moltype=211: 2'deoxy-5'tiocytidine             (r211)
%moltype=212: 2',3'-dideoxy-didehydro-cytidine           (r212)
%moltype=220: 2'deoxy-6flourinecytidine d6FC    (r220)
%moltype=230: dCMP 5'                           (r230)
%moltype=231: dCMP 3'                           (r231)
%moltype=240: rCyd                              (r240)
%moltype=250: dCyd imino                        (r250)
%moltype=270: 2'-deoxy-5-methylcytidine dm5Cyd  (r270)           2do


%moltype=301: Guanine                           (r301)
%moltype=310: 2'deoxy-8azaguanosine     d8AG    (r310)
%moltype=311: 2'deoxy-5'tioguanosine            (r311)
%moltype=312: 2',3'-dideoxy-didehydro-guanosine           (r312)
%moltype=330: dGMP 5'                           (r330)
%moltype=331: dGMP 3'                           (r331)
%moltype=340: rGuo                              (r340)
%moltype=350: dGuo lactim                       (r350)
%moltype=352: 2'deoxy-8oxo-guanosine     d8OG   (r352)
%moltype=353: 2'deoxy-7methyl-guanosine  dm7G   (r353)
%moltype=390: 2'deoxy-xanthine           dX     (r390)
%moltype=391: 2'deoxy-hypoxanthine       dH     (r391)

%moltype=401: Thymine                           (r401)
%moltype=402: Thymine crystal cell              (r402)
%moltype=411: 2'deoxy-5'tiothymidine            (r411)
%moltype=412: 2',3'-dideoxy-didehydro-thymidine           (r412) allien
%moltype=420: 2'deoxy-6flourinethymidine d6FT   (r420)
%moltype=430: dTMP 5'                           (r430)
%moltype=431: dTMP 3'                           (r431)
%moltype=450: dThd lactim??                     (r450)

%moltype=501: Uridine                           (r501)
%moltype=512: 2',3'-dideoxy-uridine           (r512)
%moltype=520: 2'deoxy-6flourineuridine  d6FU    (r520)
%moltype=521: 2'deoxy-5bromuridine  d5BrU       (r521) isakova done
%moltype=522: 2'deoxy-5bromuridine  d5BrU enol  (r522)
%moltype=530: dUMP 5'                           (r530)
%moltype=531: dUMP 3'                           (r531)
%moltype=540: rUrd                              (r540)
%moltype=550: dUrd lactim                       (r550)

%moltype=910:  double 1'2'-deoxyribose residue


%status=0: normal execution
%status~=0: error occured

clear i*

if nargin<2
    moltype=7;
end
if nargin<3
    fl_restart=0;
end


global pind
pindsdef

status=999;
lastwarn('error: identmol');

if ~fl_restart
    ms0.marked=zeros(ms0.atomnum,1,'uint16');
    ms0.pind=zeros(ms0.atomnum,1,'uint16');
else
% otherwise let leave pinds unchanged
end


if ~any(moltype==[30,101,201,301,401,501]) %not bases

    %finding O4
    iO4=0;

    if moltype>100 && mod(moltype,100)==11 %tionucleosides
      ind = findatom(ms0,'S',1);
    else
      ind = findatom(ms0,'O',1);
    end

    while ind~=0
      idxs = findbonds(ms0,ind);
      if numel(strcmpcellar(ms0.labels(idxs),'C'))==2

        if fl_restart && ms0.marked(ind)==1
           continue %let skip this atom if with are continue identification
        end

        iO4 = ind;
        if moltype>100 && mod(moltype,100)==11 %tionucleosides
          ms0.pind(ind) = find(strcmp(pind.labels,'pS4'));
        else
          ms0.pind(ind) = find(strcmp(pind.labels,'pO4'));
        end

        ms0.marked(ind)=1;
        break
      end
      ind = findatom(ms0,'O',ind+1);
    end
    if ~iO4, lastwarn('error: identmol: O4 not found'); return, end

    %finding C1
    iC1=0;
    ind = findatom(ms0,'C',1);
    while ind~=0
      idxs = findbonds(ms0,ind);
      if numel(strcmpcellar(ms0.labels(idxs),'C'))==1 && sum(idxs==iO4)==1
        iC1 = ind;
        ms0.pind(ind) = find(strcmp(pind.labels,'pC1'));
        ms0.marked(ind)=1;
        break
      end
      ind = findatom(ms0,'C',ind+1);
    end
    if ~iC1, lastwarn('error: identmol: C1 not found'); return, end

    %finding C4
%     iC4=0;
%     ind = findatom(ms0,'C',1);
%     while ind~=0
%       idxs = findbonds(ms0,ind);
%       if numel(strcmpcellar(ms0.labels(idxs),'C'))==2 && sum(idxs==iO4)==1
%         iC4 = ind;
%         ms0.pind(ind) = find(strcmp(pind.labels,'pC4'));
%         ms0.marked(ind)=1;
%         break
%       end
%       ind = findatom(ms0,'C',ind+1);
%     end
%     if ~iC4, lastwarn('error: identmol: C4 not found'); return, end
    iC4 = identatom(ms0,'C',iO4);
    ms0.pind(iC4) = find(strcmp(pind.labels,'pC4'));
    ms0.marked(iC4)=1;
    
    iC2 = identatom(ms0,'C',iC1);
    ms0.pind(iC2) = find(strcmp(pind.labels,'pC2'));
    ms0.marked(iC2)=1;

    iC3 = identatom(ms0,'C',iC2);
    if ~iC3, lastwarn('error: identmol: C3 not found IVAN');return, end %IVAN
    ms0.pind(iC3) = find(strcmp(pind.labels,'pC3'));
    ms0.marked(iC3)=1;


    if ~any(moltype==[2,3]) 
        iC5 = identatom(ms0,'C',iC4);
        if ~iC5, lastwarn('error: identmol: C5 not found IVAN');return, end %IVAN
        ms0.pind(iC5) = find(strcmp(pind.labels,'pC5'));
        ms0.marked(iC5)=1;
    end

    if any(moltype==[3,5,6,7,14]) || mod(moltype,100)==40 %ribo or 3'deoxy molecules
        iO2 = identatom(ms0,'O',iC2);
        ms0.pind(iO2) = find(strcmp(pind.labels,'pO2'));
        ms0.marked(iO2)=1;
    end

    if any(moltype~=6) && ~(moltype>100 && mod(moltype,100)==12) %3'deoxy molecules
        iO3 = identatom(ms0,'O',iC3);
        if ~iO3, lastwarn('error: identmol: O3 not found IVAN');return, end %IVAN
        ms0.pind(iO3) = find(strcmp(pind.labels,'pO3'));
        ms0.marked(iO3)=1;
    end

    if any(moltype==[4,5]) 
        %
    elseif any(moltype==[2,3]) 
        %
    else
        iO5 = identatom(ms0,'O',iC5);
        if ~iO5, lastwarn('error: identmol: O5 not found');return, end
        ms0.pind(iO5) = find(strcmp(pind.labels,'pO5'));
        ms0.marked(iO5)=1;
    end

    if any( moltype==[3,5,6,7,14]) || (moltype>100 && (mod(moltype,100)==40)) %ribo (non 2'deoxy) molecules
        iH22 = identatom(ms0,'H',iO2);
        ms0.pind(iH22) = find(strcmp(pind.labels,'pH22'));
        ms0.marked(iH22)=1;

        iH21 = identatom(ms0,'H',iC2);
        ms0.pind(iH21) = find(strcmp(pind.labels,'pH21'));
        ms0.marked(iH21)=1;

    else
%    elseif any(moltype==[8,9,11,12,13,15,16,17,18,19,21,22,23,24,25]) ||...
%       (moltype>100 && any(mod(moltype,100)==[10,20,30,31,50])) %in 1'2'deoxyribose H22 is bonded to C2

        iH21 = identatom(ms0,'H',iC2);
        if ~iH21, lastwarn('error: identmol: H21 not found IVAN');return, end %IVAN
        ms0.marked(iH21)=1;
        
        if (moltype>100 && mod(moltype,100)==12) %dideoxy molecules
            ms0.pind(iH21) = find(strcmp(pind.labels,'pH21'));
        else
            iH22 = identatom(ms0,'H',iC2);
            if ~iH22, lastwarn('error: identmol: H22 not found IVAN');return, end %IVAN
            ms0.marked(iH22)=1;

            v1 = cross(createvect(ms0,iC2,ms0,iC1),createvect(ms0,iC2,ms0,iC3));
            v2 = createvect(ms0,iC2,ms0,iH21);
            %if H21 is not upside of C1-C2-C3 plane than exchange H21 and H22 indexes
        %!!by Seanger H atom connected to O2' in ribose in 2'deoxyribose is named H21' and another ones is named H22' !! 
            if abs(p2pangle(v1,v2)) > pi/2 
              buf = iH21; iH21 = iH22; iH22 = buf;
            end

            ms0.pind(iH21) = find(strcmp(pind.labels,'pH21'));
            ms0.pind(iH22) = find(strcmp(pind.labels,'pH22'));
        end

    end

    if (moltype>100 && mod(moltype,100)==31) || moltype==23 %3'nucleotides
        ms0=detectPchain(ms0,iO3,'pP3');
    elseif (moltype>100 && mod(moltype,100)==12)  %2',3'-deoxy molecules
    elseif moltype==6  %3'deoxy molecules
        iH32 = identatom(ms0,'H',iC3);
        ms0.pind(iH32) = find(strcmp(pind.labels,'pH32'));
        ms0.marked(iH32)=1;
    else
    %!!by Seanger H atom connected to O3' in ribose in 2'deoxyribose is named H31' and another ones is named H32' !! 
        iH32 = identatom(ms0,'H',iO3);
        if ~iH32, lastwarn('error: identmol: iH32 not found IVAN');return, end %IVAN
        ms0.pind(iH32) = find(strcmp(pind.labels,'pH32'));
        ms0.marked(iH32)=1;
    end

    iH31 = identatom(ms0,'H',iC3);
    if ~iH31, lastwarn('error: identmol: iH31 not found IVAN');return, end %IVAN
    ms0.pind(iH31) = find(strcmp(pind.labels,'pH31'));
    ms0.marked(iH31)=1;

    if (moltype>100 && mod(moltype,100)==30) || moltype==22 %5'nucleotides
        ms0=detectPchain(ms0,iO5,'pP5');
    elseif any(moltype==[910])

        iP = identatom(ms0,'P',iO5);
        if ~iP, lastwarn('error: identmol: iP not found');return, end
        ms0.pind(iP) = find(strcmp(pind.labels,'pP5'));
        ms0.marked(iP)=1;
        
        iO3_2 = identatom(ms0,'O',iP,'C');
        if ~iO3_2, lastwarn('error: identmol: iO3_2 not found');return, end
        ms0.pind(iO3_2) = find(strcmp(pind.labels,'pO3'));
        ms0.marked(iO3_2)=1;
        
        iC3_2 = identatom(ms0,'C',iO3_2);
        if ~iC3_2, lastwarn('error: identmol: iC3_2 not found');return, end
        ms0.pind(iC3_2) = find(strcmp(pind.labels,'pC3'));
        ms0.marked(iC3_2)=1;
        
        iC4_2 = identatom(ms0,'C',iC3_2,'O');
        if ~iC4_2, lastwarn('error: identmol: iC4_2 not found');return, end
        ms0.pind(iC4_2) = find(strcmp(pind.labels,'pC4'));
        ms0.marked(iC4_2)=1;
        
        iO4_2 = identatom(ms0,'O',iC4_2);
        if ~iO4_2, lastwarn('error: identmol: iO4_2 not found');return, end
        ms0.pind(iO4_2) = find(strcmp(pind.labels,'pO4'));
        ms0.marked(iO4_2)=1;

        iC1_2 = identatom(ms0,'C',iO4_2);
        if ~iC1_2, lastwarn('error: identmol: iC1_2 not found');return, end
        ms0.pind(iC1_2) = find(strcmp(pind.labels,'pC1'));
        ms0.marked(iC1_2)=1;
        
        iC2_2 = identatom(ms0,'C',iC1_2);
        if ~iC2_2, lastwarn('error: identmol: iC2_2 not found');return, end
        ms0.pind(iC2_2) = find(strcmp(pind.labels,'pC2'));
        ms0.marked(iC2_2)=1;

        iC5_2 = identatom(ms0,'C',iC4_2);
        if ~iC5_2, lastwarn('error: identmol: iC5_2 not found');return, end
        ms0.pind(iC5_2) = find(strcmp(pind.labels,'pC5'));
        ms0.marked(iC5_2)=1;
        
        iO5_2 = identatom(ms0,'O',iC5_2);
        if ~iO5_2, lastwarn('error: identmol: iO5_2 not found');return, end
        ms0.pind(iO5_2) = find(strcmp(pind.labels,'pO5'));
        ms0.marked(iO5_2)=1;
        
        iH4_2 = identatom(ms0,'H',iC4_2);
        if ~iH4_2, lastwarn('error: identmol: iH4_2 not found');return, end
        ms0.pind(iH4_2) = find(strcmp(pind.labels,'pH4'));
        ms0.marked(iH4_2)=1;
        
        iH53_2 = identatom(ms0,'H',iO5_2);
        if ~iH53_2, lastwarn('error: identmol: H53_2 not found');return, end
        ms0.pind(iH53_2) = find(strcmp(pind.labels,'pH53'));
        ms0.marked(iH53_2)=1;
        
        iH31_2 = identatom(ms0,'H',iC3_2);
        if ~iH31_2, lastwarn('error: identmol: iH31_2 not found');return, end
        ms0.pind(iH31_2) = find(strcmp(pind.labels,'pH31'));
        ms0.marked(iH31_2)=1;
        
        iH51 = identatom(ms0,'H',iC5_2);
        ms0.marked(iH51)=1;
        if ~iH51, lastwarn('error: identmol: iH51 not found');return, end
        iH52 = identatom(ms0,'H',iC5_2);
        if ~iH52, lastwarn('error: identmol: iH52 not found');return, end
        ms0.marked(iH52)=1;

        v1 = cross(createvect(ms0,iC5_2,ms0,iC4_2),createvect(ms0,iC5_2,ms0,iO5_2));
        v2 = createvect(ms0,iC5_2,ms0,iH51);
        if dot(v1,v2) < 0
            buf = iH51; iH51 = iH52; iH52 = buf;
        end
        ms0.pind(iH51) = find(strcmp(pind.labels,'pH51'));
        ms0.pind(iH52) = find(strcmp(pind.labels,'pH52'));
        
        iH21 = identatom(ms0,'H',iC2_2);
        if ~iH21, lastwarn('error: identmol: iH21 not found');return, end
        ms0.marked(iH21)=1;
        iH22 = identatom(ms0,'H',iC2_2);
        if ~iH22, lastwarn('error: identmol: iH22 not found');return, end
        ms0.marked(iH22)=1;

        v1 = cross(createvect(ms0,iC2_2,ms0,iC1_2),createvect(ms0,iC2_2,ms0,iC3_2));
        v2 = createvect(ms0,iC2_2,ms0,iH21);
        %if H21 is not upside of C1-C2-C3 plane than exchange H21 and H22 indexes
        %!!by Seanger H atom connected to O2' in ribose in 2'deoxyribose is named H21' and another ones is named H22' !! 
        if abs(p2pangle(v1,v2)) > pi/2 
             buf = iH21; iH21 = iH22; iH22 = buf;
        end
        ms0.pind(iH21) = find(strcmp(pind.labels,'pH21'));
        ms0.pind(iH22) = find(strcmp(pind.labels,'pH22'));

        iH11 = identatom(ms0,'H',iC1_2);
        if ~iH11, lastwarn('error: identmol: iH11 not found');return, end
        ms0.marked(iH11)=1;
        iH12 = identatom(ms0,'H',iC1_2);
        if ~iH12, lastwarn('error: identmol: iH12 not found');return, end
        ms0.marked(iH12)=1;
        v1 = cross(createvect(ms0,iC1_2,ms0,iO4_2),createvect(ms0,iC1_2,ms0,iC2_2));
        v2 = createvect(ms0,iC1_2,ms0,iH11);
        %if H11 is not upside of O4-C1-C2 plane than exchange H11 and H12 indexes
        if abs(p2pangle(v1,v2)) > pi/2 
          buf = iH11; iH11 = iH12; iH12 = buf;
        end
        ms0.pind(iH11) = find(strcmp(pind.labels,'pH11'));
        ms0.pind(iH12) = find(strcmp(pind.labels,'pH12'));

        iOP1 = identatom(ms0,'O',iP);
        if ~iOP1, lastwarn('error: identmol: iOP1 not found');return, end
        ms0.marked(iOP1)=1;
        iOP2 = identatom(ms0,'O',iP);
        if ~iOP2, lastwarn('error: identmol: iOP2 not found');return, end
        ms0.marked(iOP2)=1;
        v1 = cross(createvect(ms0,iP,ms0,iO5),createvect(ms0,iP,ms0,iO3_2));
        v2 = createvect(ms0,iP,ms0,iOP1);
        if dot(v1,v2) < 0
            buf = iOP1; iOP1 = iOP2; iOP2 = buf;
        end
        ms0.pind(iOP1) = find(strcmp(pind.labels,'pOP1'));
        ms0.pind(iOP2) = find(strcmp(pind.labels,'pOP2'));
        
        
    elseif any(moltype==[2,3])
        % 
    elseif any(moltype==[4,5]) 
        iH53 = identatom(ms0,'H',iC5);
        ms0.pind(iH53) = find(strcmp(pind.labels,'pH53'));
        ms0.marked(iH53)=1;
    else %nucleosides
        iH53 = identatom(ms0,'H',iO5);
        if ~iH53, lastwarn('error: identmol: H53 not found');return, end
        ms0.pind(iH53) = find(strcmp(pind.labels,'pH53'));
        ms0.marked(iH53)=1;
    end

    iH4 = identatom(ms0,'H',iC4);
    ms0.pind(iH4) = find(strcmp(pind.labels,'pH4'));
    ms0.marked(iH4)=1;

    if any(moltype==[2,3])
        iH42 = identatom(ms0,'H',iC4);
        ms0.pind(iH42) = find(strcmp(pind.labels,'pH42'));
        ms0.marked(iH42)=1;
    end
        
    if ~any(moltype==[2,3])
        iH51 = identatom(ms0,'H',iC5);
        ms0.marked(iH51)=1;
        iH52 = identatom(ms0,'H',iC5);
        ms0.marked(iH52)=1;

        %!! by Seanger: when looking along O5'C5' and going clockwise we are passing atoms in following order C4' - H51' - H52' !!
    %    v1 = cross(createvect(ms0,iC5,ms0,iO5),createvect(ms0,iC5,ms0,iC4));
        if ~any(moltype==[4,5]) % r4,r5 H51 and H52 are equivalent
            v1 = cross(createvect(ms0,iC5,ms0,iC4),createvect(ms0,iC5,ms0,iO5));
            v2 = createvect(ms0,iC5,ms0,iH51);
            if dot(v1,v2) < 0
              buf = iH51; iH51 = iH52; iH52 = buf;
            end
        else
        end
        ms0.pind(iH51) = find(strcmp(pind.labels,'pH51'));
        ms0.pind(iH52) = find(strcmp(pind.labels,'pH52'));
    end

    %H11 H12
    iH12 = identatom(ms0,'H',iC1);
    ms0.marked(iH12)=1;

elseif any(moltype==30) %moltype=30: 6azacytidine base (only simple heterocyclic base)

    %finding bN1
    ibN1=0;
    ind = findatom(ms0,'N',1);
    while ind~=0
      idxs = findbonds(ms0,ind);
      if numel(strcmpcellar(ms0.labels(idxs),'N'))==1 && numel(strcmpcellar(ms0.labels(idxs)=='C'))==1
        ibN1 = ind;
        ms0.pind(ind) = find(strcmp(pind.labels,'bN1'));
        ms0.marked(ind)=1;
        break
      end
      ind = findatom(ms0,'N',ind+1);
    end
    if ~ibN1, lastwarn('error: identmol: bN1 not found'); return, end

    ibH1 = identatom(ms0,'H',ibN1);
    ms0.pind(ibH1) = find(strcmp(pind.labels,'bH1'));
    ms0.marked(ibH1)=1;

    ibC2 = identatom(ms0,'C',ibN1,'O'); % bC2 have bond with O in contrast to C6
    ms0.pind(ibC2) = find(strcmp(pind.labels,'bC2'));
    ms0.marked(ibC2)=1;

    ibN3 = identatom(ms0,'N',ibC2);
    ms0.pind(ibN3) = find(strcmp(pind.labels,'bN3'));
    ms0.marked(ibN3)=1;

    ibC4 = identatom(ms0,'C',ibN3);
    ms0.pind(ibC4) = find(strcmp(pind.labels,'bC4'));
    ms0.marked(ibC4)=1;

    ibC5 = identatom(ms0,'C',ibC4);
    ms0.pind(ibC5) = find(strcmp(pind.labels,'bC5'));
    ms0.marked(ibC5)=1;
    ibH5 = identatom(ms0,'H',ibC5);
    ms0.pind(ibH5) = find(strcmp(pind.labels,'bH5'));
    ms0.marked(ibH5)=1;

    ibN6 = identatom(ms0,'N',ibC5);
    ms0.pind(ibN6) = find(strcmp(pind.labels,'bN6'));
    ms0.marked(ibN6)=1;

    ibN = identatom(ms0,'N',ibC4);
    ms0.pind(ibN) = find(strcmp(pind.labels,'bN'));
    ms0.marked(ibN)=1;

    ibO2 = identatom(ms0,'O',ibC2);
    ms0.pind(ibO2) = find(strcmp(pind.labels,'bO2'));
    ms0.marked(ibO2)=1;

    ibH41 = identatom(ms0,'H',ibN);
    ms0.pind(ibH41) = find(strcmp(pind.labels,'bH41'));
    ms0.marked(ibH41)=1;

    ibH42 = identatom(ms0,'H',ibN);
    ms0.pind(ibH42) = find(strcmp(pind.labels,'bH42'));
    ms0.marked(ibH42)=1;

end

if any(moltype==[2,3,4,5,6,7,8,22,23,910]) %sugars
    iH11 = identatom(ms0,'H',iC1);
    ms0.marked(iH11)=1;

    v1 = cross(createvect(ms0,iC1,ms0,iO4),createvect(ms0,iC1,ms0,iC2));
    v2 = createvect(ms0,iC1,ms0,iH11);
    %if H11 is not upside of O4-C1-C2 plane than exchange H11 and H12 indexes
    if abs(p2pangle(v1,v2)) > pi/2 
      buf = iH11; iH11 = iH12; iH12 = buf;
    end

    ms0.pind(iH11) = find(strcmp(pind.labels,'pH11'));
    ms0.pind(iH12) = find(strcmp(pind.labels,'pH12'));

elseif any(moltype==15) || (moltype>=100 && moltype <= 199)%adenine

    if any(moltype==101) %adenine

      %finding bN9
      ibN9=0;
      ind = findatom(ms0,'N',1);
      while ind~=0
        idxs = findbonds(ms0,ind);
        if numel(strcmpcellar(ms0.labels(idxs),'C'))==2 && numel(strcmpcellar(ms0.labels(idxs),'H'))==1
          ibN9 = ind;
          ms0.pind(ind) = find(strcmp(pind.labels,'bN9'));
          ms0.marked(ind)=1;

          break
        end
        ind = findatom(ms0,'N',ind+1);
      end
      if ~ibN9, lastwarn('error: identmol: bN9 not found'); return, end

      ibH9 = identatom(ms0,'H',ibN9);
      ms0.pind(ibH9) = find(strcmp(pind.labels,'bH9'));
      ms0.marked(ibH9)=1;

    else

      ms0.pind(iH12) = find(strcmp(pind.labels,'pH12'));

      ibN9 = identatom(ms0,'N',iC1);
      ms0.pind(ibN9) = find(strcmp(pind.labels,'bN9'));
      ms0.marked(ibN9)=1;
    end

    ibC4 = identatom(ms0,'C',ibN9,'C'); % bC4 have bond with C5 in contrast to C8
    ms0.pind(ibC4) = find(strcmp(pind.labels,'bC4'));
    ms0.marked(ibC4)=1;

    ibN3 = identatom(ms0,'N',ibC4);
    ms0.pind(ibN3) = find(strcmp(pind.labels,'bN3'));
    ms0.marked(ibN3)=1;

    ibC2 = identatom(ms0,'C',ibN3);
    ms0.pind(ibC2) = find(strcmp(pind.labels,'bC2'));
    ms0.marked(ibC2)=1;


    if moltype==151 %d2aminopurine

        ibCN2 = identatom(ms0,'N',ibC2,'H');
        ms0.pind(ibCN2) = find(strcmp(pind.labels,'bCN2'));
        ms0.marked(ibCN2)=1;
        ibH21 = identatom(ms0,'H',ibCN2);
        ms0.pind(ibH21) = find(strcmp(pind.labels,'bH21'));
        ms0.marked(ibH21)=1;

        ibH22 = identatom(ms0,'H',ibCN2);
        ms0.pind(ibH22) = find(strcmp(pind.labels,'bH22'));
        ms0.marked(ibH22)=1;
    else
        ibH2 = identatom(ms0,'H',ibC2);
        ms0.pind(ibH2) = find(strcmp(pind.labels,'bH2'));
        ms0.marked(ibH2)=1;
    end

    ibN1 = identatom(ms0,'N',ibC2);
    ms0.pind(ibN1) = find(strcmp(pind.labels,'bN1'));
    ms0.marked(ibN1)=1;

    ibC6 = identatom(ms0,'C',ibN1);
    ms0.pind(ibC6) = find(strcmp(pind.labels,'bC6'));
    ms0.marked(ibC6)=1;

    if any(moltype==[151,161,162]) %d2aminopurine, 2'-дезоксипурин, 2'-дезокси-8-aza-пурин
        ibH6 = identatom(ms0,'H',ibC6);
        ms0.pind(ibH6) = find(strcmp(pind.labels,'bH6'));
        ms0.marked(ibH6)=1;
    else

        ibCN6 = identatom(ms0,'N',ibC6);
        ms0.pind(ibCN6) = find(strcmp(pind.labels,'bCN6'));
        ms0.marked(ibCN6)=1;
        ibH61 = identatom(ms0,'H',ibCN6);
        ms0.pind(ibH61) = find(strcmp(pind.labels,'bH61'));
        ms0.marked(ibH61)=1;

        if moltype==150 %dAde imino
            ibH1 = identatom(ms0,'H',ibN1);
            ms0.pind(ibH1) = find(strcmp(pind.labels,'bH1'));
            ms0.marked(ibH1)=1;
        else
            ibH62 = identatom(ms0,'H',ibCN6);
            ms0.pind(ibH62) = find(strcmp(pind.labels,'bH62'));
            ms0.marked(ibH62)=1;
        end
    end 

    ibC5 = identatom(ms0,'C',ibC6);
    ms0.pind(ibC5) = find(strcmp(pind.labels,'bC5'));
    ms0.marked(ibC5)=1;

    ibN7 = identatom(ms0,'N',ibC5);
    ms0.pind(ibN7) = find(strcmp(pind.labels,'bN7'));
    ms0.marked(ibN7)=1;

    if any(moltype==[110,162]) %8azaadenosine, 2'-дезокси-8-aza-пурин
        ibN8 = identatom(ms0,'N',ibN7);
        ms0.pind(ibN8) = find(strcmp(pind.labels,'bN8'));
        ms0.marked(ibN8)=1;
    else
        ibC8 = identatom(ms0,'C',ibN7);
        ms0.pind(ibC8) = find(strcmp(pind.labels,'bC8'));
        ms0.marked(ibC8)=1;
          
        if moltype==152 %8oxy
            ibCO8 = identatom(ms0,'O',ibC8);
            ms0.pind(ibCO8) = find(strcmp(pind.labels,'bCO8'));
            ms0.marked(ibCO8)=1;

            ibH7 = identatom(ms0,'H',ibN7);
            ms0.pind(ibH7) = find(strcmp(pind.labels,'bH7'));
            ms0.marked(ibH7)=1;
        else
            ibH8 = identatom(ms0,'H',ibC8);
            ms0.pind(ibH8) = find(strcmp(pind.labels,'bH8'));
            ms0.marked(ibH8)=1;
        end

    end

elseif any(moltype==[11,9,14,21]) || (moltype>=200 && moltype <= 299) %cytidine & azocytidine & flourinecytidine molecules

    if any(moltype==201) %cytidine

      %finding bN1
      ibN1=0;
      ind = findatom(ms0,'N',1);
      while ind~=0
        idxs = findbonds(ms0,ind);
        if numel(strcmpcellar(ms0.labels(idxs),'C'))==2 
          ibN1 = ind;

          ibC2 = identatom(ms0,'C',ibN1,'O'); % bC2 have bond with O in contrast to C6
          if ~ibC2, continue, end;

          break
        end
        ind = findatom(ms0,'N',ind+1);
      end
      if ~ibN1, lastwarn('error: identmol: bN1 not found'); return, end

      ms0.pind(ind) = find(strcmp(pind.labels,'bN1'));
      ms0.marked(ind)=1;
      ms0.pind(ibC2) = find(strcmp(pind.labels,'bC2'));
      ms0.marked(ibC2)=1;

      ibH1 = identatom(ms0,'H',ibN1);
      ms0.pind(ibH1) = find(strcmp(pind.labels,'bH1'));
      ms0.marked(ibH1)=1;

    else
      ms0.pind(iH12) = find(strcmp(pind.labels,'pH12'));

      ibN1 = identatom(ms0,'N',iC1);
      ms0.pind(ibN1) = find(strcmp(pind.labels,'bN1'));
      ms0.marked(ibN1)=1;

      ibC2 = identatom(ms0,'C',ibN1,'O'); % bC2 have bond with O in contrast to C6
      ms0.pind(ibC2) = find(strcmp(pind.labels,'bC2'));
      ms0.marked(ibC2)=1;
    end

    ibN3 = identatom(ms0,'N',ibC2);
    ms0.pind(ibN3) = find(strcmp(pind.labels,'bN3'));
    ms0.marked(ibN3)=1;

    ibC4 = identatom(ms0,'C',ibN3);
    ms0.pind(ibC4) = find(strcmp(pind.labels,'bC4'));
    ms0.marked(ibC4)=1;

    if moltype==21 %5AC

      ibN5 = identatom(ms0,'N',ibC4,'C','H');
      ms0.pind(ibN5) = find(strcmp(pind.labels,'bN5'));
      ms0.marked(ibN5)=1;

      ibC6 = identatom(ms0,'C',ibN5);
      ms0.pind(ibC6) = find(strcmp(pind.labels,'bC6'));
      ms0.marked(ibC6)=1;
      ibH6 = identatom(ms0,'H',ibC6);
      ms0.pind(ibH6) = find(strcmp(pind.labels,'bH6'));
      ms0.marked(ibH6)=1;

    else %not 5AC molecules

      ibC5 = identatom(ms0,'C',ibC4);
      ms0.pind(ibC5) = find(strcmp(pind.labels,'bC5'));
      ms0.marked(ibC5)=1;

      if any(moltype==[11,14]) %d6AC or 6AC
        ibN6 = identatom(ms0,'N',ibC5);
        ms0.pind(ibN6) = find(strcmp(pind.labels,'bN6'));
        ms0.marked(ibN6)=1;
      elseif moltype==220 % 6Fluorine molecules 
        ibC6 = identatom(ms0,'C',ibC5);
        ms0.pind(ibC6) = find(strcmp(pind.labels,'bC6'));
        ms0.marked(ibC6)=1;
        ibF6 = identatom(ms0,'F',ibC6);
        ms0.pind(ibF6) = find(strcmp(pind.labels,'bF6'));
        ms0.marked(ibF6)=1;
      else %if any(moltype==[9,230,231,240,250]) %molecules that have C6H6 bond
        ibC6 = identatom(ms0,'C',ibC5);
        ms0.pind(ibC6) = find(strcmp(pind.labels,'bC6'));
        ms0.marked(ibC6)=1;
        ibH6 = identatom(ms0,'H',ibC6);
        ms0.pind(ibH6) = find(strcmp(pind.labels,'bH6'));
        ms0.marked(ibH6)=1;
      end
      
      if moltype==270 %d5MCyd
          
          [ibCC5,ibH51,ibH52,ibH53]=identCH3(ibC5,ibC6);
%          ibCC5 = identatom(ms0,'C',ibC5,'','N');
%
%          %CH3 group atoms are now distinguished by C4C5CC5Hx torsions
%          ibH51 = identatom(ms0,'H',ibCC5);
%          ibH52 = identatom(ms0,'H',ibCC5);
%          ibH53 = identatom(ms0,'H',ibCC5);
%
%          atom6 = ibC6; %index of atom at 6 position of base
%          tor(1) = torang(ms0,atom6,ibC5,ibCC5,ibH51);
%          tor(2) = torang(ms0,atom6,ibC5,ibCC5,ibH52);
%          tor(3) = torang(ms0,atom6,ibC5,ibCC5,ibH53);
%          ibH5(1) = ibH51;  ibH5(2) = ibH52;  ibH5(3) = ibH53;
%
%          [mintor,i] = min(abs(tor));
%          if mintor>30
%            warning('identmol:dThd','CH3 dThd error: C4C5CC5H1 torsion is larger than 30 degree. Improper CH3 location');
%          end
%          ibH51=ibH5(i);
%          ibH52 = setdiff(ibH5(tor>0),ibH5(i));
%          ibH53 = setdiff(ibH5(tor<0),ibH5(i));

          ms0.pind(ibCC5) = find(strcmp(pind.labels,'bCC5'));
%          ms0.marked(ibCC5)=1;

 %         ms0.marked(ibH51)=1;
 %         ms0.marked(ibH52)=1;
 %         ms0.marked(ibH53)=1;
          ms0.pind(ibH51) = find(strcmp(pind.labels,'bH51')); 
          ms0.pind(ibH52) = find(strcmp(pind.labels,'bH52'));
          ms0.pind(ibH53) = find(strcmp(pind.labels,'bH53'));
      else %not d5MCyd
          ibH5 = identatom(ms0,'H',ibC5);
          ms0.pind(ibH5) = find(strcmp(pind.labels,'bH5'));
          ms0.marked(ibH5)=1;
      end

    end


    ibN = identatom(ms0,'N',ibC4);
    ms0.pind(ibN) = find(strcmp(pind.labels,'bN'));
    ms0.marked(ibN)=1;

    ibO2 = identatom(ms0,'O',ibC2);
    ms0.pind(ibO2) = find(strcmp(pind.labels,'bO2'));
    ms0.marked(ibO2)=1;

    ibH41 = identatom(ms0,'H',ibN);
    ms0.pind(ibH41) = find(strcmp(pind.labels,'bH41'));
    ms0.marked(ibH41)=1;

    if moltype==250 %dCyd imino
        ibH3 = identatom(ms0,'H',ibN3);
        ms0.pind(ibH3) = find(strcmp(pind.labels,'bH3'));
        ms0.marked(ibH3)=1;
    else
        ibH42 = identatom(ms0,'H',ibN);
        ms0.pind(ibH42) = find(strcmp(pind.labels,'bH42'));
        ms0.marked(ibH42)=1;
    end
    %cytidine based molecules analysis end

elseif any(moltype==16) || (moltype>=300 && moltype <= 399)%guanine %ksantine

    if any(moltype==301) %guanine

      %finding bN9
      ibN9=0;
      ind = findatom(ms0,'N',1);
      while ind~=0
        idxs = findbonds(ms0,ind);
        if numel(strcmpcellar(ms0.labels(idxs),'C'))==2 && numel(strcmpcellar(ms0.labels(idxs),'H'))==1
          ibN9 = ind;

          %ixC = identatom(ms0,'C',ibN9);
          if numel(strcmpcellar(ms0.labels(idxs),'N'))==3 %test if this is not bC2 atom
            continue
          end
          if numel(strcmpcellar(ms0.labels(idxs),'O'))==1 %test if this is not bC6 atom
            continue
          end
          ms0.pind(ind) = find(strcmp(pind.labels,'bN9'));
          ms0.marked(ind)=1;

          break
        end
        ind = findatom(ms0,'N',ind+1);
      end
      if ~ibN9, lastwarn('error: identmol: bN9 not found'); return, end

      ibH9 = identatom(ms0,'H',ibN9);
      ms0.pind(ibH9) = find(strcmp(pind.labels,'bH9'));
      ms0.marked(ibH9)=1;

    else

      ms0.pind(iH12) = find(strcmp(pind.labels,'pH12'));

      ibN9 = identatom(ms0,'N',iC1);
      ms0.pind(ibN9) = find(strcmp(pind.labels,'bN9'));
      ms0.marked(ibN9)=1;
    end

    ibC4 = identatom(ms0,'C',ibN9,'C'); % bC4 have bond with C5 in contrast to C8
    ms0.pind(ibC4) = find(strcmp(pind.labels,'bC4'));
    ms0.marked(ibC4)=1;

    ibN3 = identatom(ms0,'N',ibC4);
    ms0.pind(ibN3) = find(strcmp(pind.labels,'bN3'));
    ms0.marked(ibN3)=1;

    if any(moltype==390) %hipoksantine
        ibH3 = identatom(ms0,'H',ibN3);
        ms0.pind(ibH3) = find(strcmp(pind.labels,'bH3'));
        ms0.marked(ibH3)=1;
    end

    ibC2 = identatom(ms0,'C',ibN3);
    ms0.pind(ibC2) = find(strcmp(pind.labels,'bC2'));
    ms0.marked(ibC2)=1;

    ibN1 = identatom(ms0,'N',ibC2,'C'); % bN1 have bond with C6 in contrast to CN2
    ms0.pind(ibN1) = find(strcmp(pind.labels,'bN1'));
    ms0.marked(ibN1)=1;

    if any(moltype==390) %ksantine
        ibO2 = identatom(ms0,'O',ibC2);
        ms0.pind(ibO2) = find(strcmp(pind.labels,'bO2'));
        ms0.marked(ibO2)=1;
    elseif any(moltype==391) %hipoksantine
        ibH2 = identatom(ms0,'H',ibC2);
        ms0.pind(ibH2) = find(strcmp(pind.labels,'bH2'));
        ms0.marked(ibH2)=1;
    else
        ibCN2 = identatom(ms0,'N',ibC2);
        ms0.pind(ibCN2) = find(strcmp(pind.labels,'bCN2'));
        ms0.marked(ibCN2)=1;
        ibH21 = identatom(ms0,'H',ibCN2);
        ms0.pind(ibH21) = find(strcmp(pind.labels,'bH21'));
        ms0.marked(ibH21)=1;
        ibH22 = identatom(ms0,'H',ibCN2);
        ms0.pind(ibH22) = find(strcmp(pind.labels,'bH22'));
        ms0.marked(ibH22)=1;
    end

    ibC6 = identatom(ms0,'C',ibN1);
    ms0.pind(ibC6) = find(strcmp(pind.labels,'bC6'));
    ms0.marked(ibC6)=1;
    ibCO6 = identatom(ms0,'O',ibC6);
    ms0.pind(ibCO6) = find(strcmp(pind.labels,'bCO6'));
    ms0.marked(ibCO6)=1;

    ibC5 = identatom(ms0,'C',ibC6);
    ms0.pind(ibC5) = find(strcmp(pind.labels,'bC5'));
    ms0.marked(ibC5)=1;

    ibN7 = identatom(ms0,'N',ibC5);
    ms0.pind(ibN7) = find(strcmp(pind.labels,'bN7'));
    ms0.marked(ibN7)=1;

    if moltype==310 %8aza
        ibN8 = identatom(ms0,'N',ibN7);
        ms0.pind(ibN8) = find(strcmp(pind.labels,'bN8'));
        ms0.marked(ibN8)=1;
    else
        ibC8 = identatom(ms0,'C',ibN7);
        ms0.pind(ibC8) = find(strcmp(pind.labels,'bC8'));
        ms0.marked(ibC8)=1;
          
        if moltype==352 %8oxy
            ibCO8 = identatom(ms0,'O',ibC8);
            ms0.pind(ibCO8) = find(strcmp(pind.labels,'bCO8'));
            ms0.marked(ibCO8)=1;

            ibH7 = identatom(ms0,'H',ibN7);
            ms0.pind(ibH7) = find(strcmp(pind.labels,'bH7'));
            ms0.marked(ibH7)=1;
        else
            ibH8 = identatom(ms0,'H',ibC8);
            ms0.pind(ibH8) = find(strcmp(pind.labels,'bH8'));
            ms0.marked(ibH8)=1;
        end

    end

    if any(moltype==353) %7-methyl-guanosine
        [ibNC7,ibH1,ibH2,ibH3]=identCH3(ibN7,ibC8);
        ms0.pind(ibNC7) = find(strcmp(pind.labels,'bNC7'));
%        ms0.marked(ibNC7)=1;

 %       ms0.marked(ibH1)=1;
 %       ms0.marked(ibH2)=1;
 %       ms0.marked(ibH3)=1;
        ms0.pind(ibH1) = find(strcmp(pind.labels,'bH71')); 
        ms0.pind(ibH2) = find(strcmp(pind.labels,'bH72'));
        ms0.pind(ibH3) = find(strcmp(pind.labels,'bH73'));
    end

    if moltype==350 %dGuo lactim
        ibH6 = identatom(ms0,'H',ibCO6);
        ms0.pind(ibH6) = find(strcmp(pind.labels,'bH6'));
        ms0.marked(ibH6)=1;
    end

    if all(moltype~=[350,353]) %dGuo lactim, 7-methyl
        ibH1 = identatom(ms0,'H',ibN1);
        ms0.pind(ibH1) = find(strcmp(pind.labels,'bH1'));
        ms0.marked(ibH1)=1;
    end

elseif any(moltype==[13,18]) || (moltype>=400 && moltype <= 499)%thymidine

    if any(moltype==401) %thymine

      %finding bN1
      ibN1=0;
      ind = findatom(ms0,'N',1);
      while ind~=0
        idxs = findbonds(ms0,ind);
        if numel(strcmpcellar(ms0.labels(idxs),'C'))==2 
          ibN1 = ind;

          ibC2 = identatom(ms0,'C',ibN1,'O'); % bC2 have bond with O in contrast to C6
          if ~ibC2, continue, end;

          break
        end
        ind = findatom(ms0,'N',ind+1);
      end
      if ~ibN1, lastwarn('error: identmol: bN1 not found'); return, end

      ms0.pind(ind) = find(strcmp(pind.labels,'bN1'));
      ms0.marked(ind)=1;
      ms0.pind(ibC2) = find(strcmp(pind.labels,'bC2'));
      ms0.marked(ibC2)=1;

      ibH1 = identatom(ms0,'H',ibN1);
      ms0.pind(ibH1) = find(strcmp(pind.labels,'bH1'));
      ms0.marked(ibH1)=1;

    else
      ms0.pind(iH12) = find(strcmp(pind.labels,'pH12'));

      ibN1 = identatom(ms0,'N',iC1);
      ms0.pind(ibN1) = find(strcmp(pind.labels,'bN1'));
      ms0.marked(ibN1)=1;

      ibC2 = identatom(ms0,'C',ibN1,'O'); % bC2 have bond with O in contrast to C6
      ms0.pind(ibC2) = find(strcmp(pind.labels,'bC2'));
      ms0.marked(ibC2)=1;
    end



    ibN3 = identatom(ms0,'N',ibC2);
    ms0.pind(ibN3) = find(strcmp(pind.labels,'bN3'));
    ms0.marked(ibN3)=1;

    ibC4 = identatom(ms0,'C',ibN3);
    ms0.pind(ibC4) = find(strcmp(pind.labels,'bC4'));
    ms0.marked(ibC4)=1;

    ibC5 = identatom(ms0,'C',ibC4);
    ms0.pind(ibC5) = find(strcmp(pind.labels,'bC5'));
    ms0.marked(ibC5)=1;

    if moltype==18 %6aza
        ibN6 = identatom(ms0,'N',ibC5);
        ms0.pind(ibN6) = find(strcmp(pind.labels,'bN6'));
        ms0.marked(ibN6)=1;
    else
        ibC6 = identatom(ms0,'C',ibN1);
        ms0.pind(ibC6) = find(strcmp(pind.labels,'bC6'));
        ms0.marked(ibC6)=1;

        if moltype==420 %6flourine
          ibF6 = identatom(ms0,'F',ibC6);
          ms0.pind(ibF6) = find(strcmp(pind.labels,'bF6'));
          ms0.marked(ibF6)=1;
        else
          ibH6 = identatom(ms0,'H',ibC6);
          ms0.pind(ibH6) = find(strcmp(pind.labels,'bH6'));
          ms0.marked(ibH6)=1;
        end
    end

    if moltype==18
        [ibCC5,ibH51,ibH52,ibH53]=identCH3(ibC5,ibN6);
    else
        [ibCC5,ibH51,ibH52,ibH53]=identCH3(ibC5,ibC6);
    end          
%    ibCC5 = identatom(ms0,'C',ibC5);
%    %CH3 group atoms are now distinguished by C4C5CC5Hx torsions
%    ibH51 = identatom(ms0,'H',ibCC5);
%    ibH52 = identatom(ms0,'H',ibCC5);
%    ibH53 = identatom(ms0,'H',ibCC5);
%    if moltype==18
%        atom6 = ibN6;
%    else
%        atom6 = ibC6; %index of atom at 6 position of base
%    end
%    tor(1) = torang(ms0,atom6,ibC5,ibCC5,ibH51);
%    tor(2) = torang(ms0,atom6,ibC5,ibCC5,ibH52);
%    tor(3) = torang(ms0,atom6,ibC5,ibCC5,ibH53);
%    ibH5(1) = ibH51;  ibH5(2) = ibH52;  ibH5(3) = ibH53;
%    [mintor,i] = min(abs(tor));
%    if mintor>30
%        warning('identmol:dThd','CH3 dThd error: C4C5CC5H1 torsion is larger than 30 degree. Improper CH3 location');
%    end
%    ibH51=ibH5(i);
%    ibH52 = setdiff(ibH5(tor>0),ibH5(i));
%    ibH53 = setdiff(ibH5(tor<0),ibH5(i));
    
    ms0.pind(ibCC5) = find(strcmp(pind.labels,'bCC5'));
%    ms0.marked(ibCC5)=1;

%    ms0.marked(ibH51)=1;
%    ms0.marked(ibH52)=1;
%    ms0.marked(ibH53)=1;
    ms0.pind(ibH51) = find(strcmp(pind.labels,'bH51')); 
    ms0.pind(ibH52) = find(strcmp(pind.labels,'bH52'));
    ms0.pind(ibH53) = find(strcmp(pind.labels,'bH53'));

    
    ibO4 = identatom(ms0,'O',ibC4);
    ms0.pind(ibO4) = find(strcmp(pind.labels,'bO4'));
    ms0.marked(ibO4)=1;

    ibO2 = identatom(ms0,'O',ibC2);
    ms0.pind(ibO2) = find(strcmp(pind.labels,'bO2'));
    ms0.marked(ibO2)=1;

    if moltype==450 %dUrd lactim
        ibH4 = identatom(ms0,'H',ibO4);
        ms0.pind(ibH4) = find(strcmp(pind.labels,'bH4'));
        ms0.marked(ibH4)=1;
    else
        ibH3 = identatom(ms0,'H',ibN3);
        ms0.pind(ibH3) = find(strcmp(pind.labels,'bH3'));
        ms0.marked(ibH3)=1;
    end

elseif any(moltype==[12,17,19]) || (moltype>=500 && moltype <= 599)%uridine

    if any(moltype==501) %uridine

      %finding bN1
      ibN1=0;
      ind = findatom(ms0,'N',1);
      while ind~=0
        idxs = findbonds(ms0,ind);
        if numel(strcmpcellar(ms0.labels(idxs),'C'))==2 
          ibN1 = ind;

          ibC2 = identatom(ms0,'C',ibN1,'O'); % bC2 have bond with O in contrast to C6
          if ~ibC2, continue, end;

          break
        end
        ind = findatom(ms0,'N',ind+1);
      end
      if ~ibN1, lastwarn('error: identmol: bN1 not found'); return, end

      ms0.pind(ind) = find(strcmp(pind.labels,'bN1'));
      ms0.marked(ind)=1;
      ms0.pind(ibC2) = find(strcmp(pind.labels,'bC2'));
      ms0.marked(ibC2)=1;

      ibH1 = identatom(ms0,'H',ibN1);
      ms0.pind(ibH1) = find(strcmp(pind.labels,'bH1'));
      ms0.marked(ibH1)=1;

    else
      ms0.pind(iH12) = find(strcmp(pind.labels,'pH12'));

      ibN1 = identatom(ms0,'N',iC1);
      ms0.pind(ibN1) = find(strcmp(pind.labels,'bN1'));
      ms0.marked(ibN1)=1;

      ibC2 = identatom(ms0,'C',ibN1,'O'); % bC2 have bond with O in contrast to C6
      ms0.pind(ibC2) = find(strcmp(pind.labels,'bC2'));
      ms0.marked(ibC2)=1;
    end

    ibN3 = identatom(ms0,'N',ibC2);
    ms0.pind(ibN3) = find(strcmp(pind.labels,'bN3'));
    ms0.marked(ibN3)=1;

    ibC4 = identatom(ms0,'C',ibN3);
    ms0.pind(ibC4) = find(strcmp(pind.labels,'bC4'));
    ms0.marked(ibC4)=1;

    if moltype==19 %5aza
        ibN5 = identatom(ms0,'N',ibC4);
        ms0.pind(ibN5) = find(strcmp(pind.labels,'bN5'));
        ms0.marked(ibN5)=1;

        ibC6 = identatom(ms0,'C',ibN5);
        ms0.pind(ibC6) = find(strcmp(pind.labels,'bC6'));
        ms0.marked(ibC6)=1;
        ibH6 = identatom(ms0,'H',ibC6);
        ms0.pind(ibH6) = find(strcmp(pind.labels,'bH6'));
        ms0.marked(ibH6)=1;
    else
        ibC5 = identatom(ms0,'C',ibC4);
        ms0.pind(ibC5) = find(strcmp(pind.labels,'bC5'));
        ms0.marked(ibC5)=1;

        if any(moltype==[521,522]) %5Br Urd
            ibBr5 = identatom(ms0,'Br',ibC5);
            ms0.pind(ibBr5) = find(strcmp(pind.labels,'bBr5'));
            ms0.marked(ibBr5)=1;
        else
            ibH5 = identatom(ms0,'H',ibC5);
            ms0.pind(ibH5) = find(strcmp(pind.labels,'bH5'));
            ms0.marked(ibH5)=1;
        end

        if moltype==17 %6aza
            ibN6 = identatom(ms0,'N',ibC5);
            ms0.pind(ibN6) = find(strcmp(pind.labels,'bN6'));
            ms0.marked(ibN6)=1;
        else
            ibC6 = identatom(ms0,'C',ibC5);
            ms0.pind(ibC6) = find(strcmp(pind.labels,'bC6'));
            ms0.marked(ibC6)=1;

            if moltype==520 %6flourine
                ibF6 = identatom(ms0,'F',ibC6);
                ms0.pind(ibF6) = find(strcmp(pind.labels,'bF6'));
                ms0.marked(ibF6)=1;
            else
                ibH6 = identatom(ms0,'H',ibC6);
                ms0.pind(ibH6) = find(strcmp(pind.labels,'bH6'));
                ms0.marked(ibH6)=1;
            end

        end
    end


    ibO4 = identatom(ms0,'O',ibC4);
    ms0.pind(ibO4) = find(strcmp(pind.labels,'bO4'));
    ms0.marked(ibO4)=1;

    if any(moltype==[550,522]) %dUrd lactim
        ibH4 = identatom(ms0,'H',ibO4);
        ms0.pind(ibH4) = find(strcmp(pind.labels,'bH4'));
        ms0.marked(ibH4)=1;
    else
        ibH3 = identatom(ms0,'H',ibN3);
        ms0.pind(ibH3) = find(strcmp(pind.labels,'bH3'));
        ms0.marked(ibH3)=1;
    end

    ibO2 = identatom(ms0,'O',ibC2);
    ms0.pind(ibO2) = find(strcmp(pind.labels,'bO2'));
    ms0.marked(ibO2)=1;

end

inds=find(ms0.pind==0);

if ~isempty(inds)
    ms0.pind(inds) = numel(pind.labels);
    status=1;
    lastwarn('Not all atoms are identified');
else
    status=0;
    lastwarn('');
end

if 0
%add information about molecule orientation in relation to internal coordinates 
%cordinate systen orientation is determined by 3 atoms oriantation: 
%Ox axis is placed at order(2)-order(1) direction
%OxOy plane is coincided with order(1)-order(2)-order(3) plane
%Oz direction is determined by Ox*Oy cross product
if any(moltype == [7 8 9 11]) > 0
    order=[6 1 2]; %'pO4' 'pC1' 'pC2'
    ind(1)=find(ms0.pind==order(1));
    ind(2)=find(ms0.pind==order(2));
    ind(3)=find(ms0.pind==order(3));
    vx=[ms0.x(ind(1))-ms0.x(ind(2)) ms0.y(ind(1))-ms0.y(ind(2)) ms0.z(ind(1))-ms0.z(ind(2))];
    ox=vx/norm(vx);
    vy=[ms0.x(ind(3))-ms0.x(ind(2)) ms0.y(ind(3))-ms0.y(ind(2)) ms0.z(ind(3))-ms0.z(ind(2))];
    vz=cross(vx,vy);
    oz=vz/norm(vz);
    oy=cross(oz,ox);
    ms0.intort=[ox;oy;oz];
ms0.intort  

else
  warning('identmol:orient','can"t to calculate intort !');
end
end

  % -----------------------------------------------------------------------
  % Nested function
  %

  function ms0=detectPchain(ms0,iOatom,Patomstring)
  %detect P residue chain connected to iOatom

    global CH
    iP = identatom(ms0,'P',iOatom);
    ms0.pind(iP) = find(strcmp(pind.labels,Patomstring));
    ms0.marked(iP)=1;

    iOPx = identatom(ms0,'O',iP,'H');
    ms0.marked(iOPx)=1;
    iHPx = identatom(ms0,'H',iOPx);
    ms0.marked(iHPx)=1;
    iOPy = identatom(ms0,'O',iP,'H');
    ms0.marked(iOPy)=1;
    iHPy = identatom(ms0,'H',iOPy);
    ms0.marked(iHPy)=1;

    [ind,dist]=findnearestatom(ms0,'N',iHPx);
    if (dist < CH.typOHNdist) && ms0.pind(ind)==0
      ixN=ind;
      fl=1;
    end
    [ind,dist]=findnearestatom(ms0,'N',iHPy);
    if (dist < CH.typOHNdist) && ms0.pind(ind)==0
      ixN=ind;
      fl=2;
    end

    if ~ixN, lastwarn('error: identmol: nucleotide analysis specified but xN atom is not found');return, end

    ms0.pind(ixN) = find(strcmp(pind.labels,'xN'));
    ms0.marked(ixN)=1;
    ms0.btA(end+1)=iP; %adding connection between nucleotide and NH3 molecule
    ms0.btB(end+1)=ixN;
    
    if fl==1  % making i*P1 atoms to be nearest to HN3 residue ones
      iOP1=iOPx; iHP1=iHPx; iOP2=iOPy; iHP2=iHPy;
    else
      iOP1=iOPy; iHP1=iHPy; iOP2=iOPx; iHP2=iHPx;
    end

    ms0.pind(iOP1) = find(strcmp(pind.labels,'pOP1'));
    ms0.pind(iHP1) = find(strcmp(pind.labels,'pHP1'));
    ms0.pind(iOP2) = find(strcmp(pind.labels,'pOP2'));
    ms0.pind(iHP2) = find(strcmp(pind.labels,'pHP2'));

    iOP3 = identatom(ms0,'O',iP);
    ms0.pind(iOP3) = find(strcmp(pind.labels,'pOP3'));
    ms0.marked(iOP3)=1;


    ixH1 = identatom(ms0,'H',ixN);
    ms0.pind(ixH1) = find(strcmp(pind.labels,'xH1'));
    ms0.marked(ixH1)=1;
    ixH2 = identatom(ms0,'H',ixN);
    ms0.pind(ixH2) = find(strcmp(pind.labels,'xH2'));
    ms0.marked(ixH2)=1;
    ixH3 = identatom(ms0,'H',ixN);
    ms0.pind(ixH3) = find(strcmp(pind.labels,'xH3'));
    ms0.marked(ixH3)=1;
  end

  function [iC,iH1,iH2,iH3] = identCH3(iA, iX)
  %identify CH3 group

    iC = identatom(ms0,'C',iA);
    ms0.marked(iC)=1;

    %CH3 group atoms are now distinguished by X-A-C-Hx torsions
    iH1 = identatom(ms0,'H',iC);
    ms0.marked(iH1)=1;
    iH2 = identatom(ms0,'H',iC);
    ms0.marked(iH2)=1;
    iH3 = identatom(ms0,'H',iC);
    ms0.marked(iH3)=1;

    tor(1) = torang(ms0,iX,iA,iC,iH1);
    tor(2) = torang(ms0,iX,iA,iC,iH2);
    tor(3) = torang(ms0,iX,iA,iC,iH3);
    iH(1) = iH1;  iH(2) = iH2;  iH(3) = iH3;

    [mintor,i] = min(abs(tor));
    if mintor>30
        [mintor,i] = min(abs(abs(tor)-180));
    end
    if mintor>30
        warning('identmol:CH3error','CH3 error: X-A-C-H1 torsion is larger than 30 degree. Improper CH3 location');
    end
    iH1 = iH(i);
    iH2 = setdiff(iH(tor>0),iH(i));
    iH3 = setdiff(iH(tor<0),iH(i));
    
  end
  % -----------------------------------------------------------------------

end %identmol
