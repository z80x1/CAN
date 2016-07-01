function ms0=calcproperties910(ms0,moltype)
%calculate molecule r910 properties
%
% Version 1.0    
% Last modified  gator 2009-08-19
% Created        gator 2009-07-?


if nargin<2
    error('moltype is not specified!');
end

global GLsdesc GLsdescr
pindsdef
atomsind
    
%ident_atom
    iC1  = ms0.ind(find(find(strcmp(pind.labels,'pC1'))==ms0.pind)); 
    iC2  = ms0.ind(find(find(strcmp(pind.labels,'pC2'))==ms0.pind)); 
    iC3  = ms0.ind(find(find(strcmp(pind.labels,'pC3'))==ms0.pind)); 
    iC4  = ms0.ind(find(find(strcmp(pind.labels,'pC4'))==ms0.pind)); 
    iC5  = ms0.ind(find(find(strcmp(pind.labels,'pC5'))==ms0.pind)); 
    iP5  = ms0.ind(find(find(strcmp(pind.labels,'pP5'))==ms0.pind));
    iO3  = ms0.ind(find(find(strcmp(pind.labels,'pO3'))==ms0.pind));        
    iO4  = ms0.ind(find(find(strcmp(pind.labels,'pO4'))==ms0.pind));    
    iO5  = ms0.ind(find(find(strcmp(pind.labels,'pO5'))==ms0.pind)); 
    iH32 = ms0.ind(find(find(strcmp(pind.labels,'pH32'))==ms0.pind)); 
    iH53 = ms0.ind(find(find(strcmp(pind.labels,'pH53'))==ms0.pind)); 
   
%sugar 1
      prop.tnu0=torang(ms0,iC4(1),iO4(1),iC1(1),iC2(1));
      prop.tnu1=torang(ms0,iO4(1),iC1(1),iC2(1),iC3(1));
      prop.tnu2=torang(ms0,iC1(1),iC2(1),iC3(1),iC4(1));
      prop.tnu3=torang(ms0,iC2(1),iC3(1),iC4(1),iO4(1));
      prop.tnu4=torang(ms0,iC3(1),iC4(1),iO4(1),iC1(1));
     prop.rnu0=adist(ms0,iO4(1),iC1(1));
     prop.rnu1=adist(ms0,iC1(1),iC2(1));
     prop.rnu2=adist(ms0,iC2(1),iC3(1));
     prop.rnu3=adist(ms0,iC3(1),iC4(1));
     prop.rnu4=adist(ms0,iC4(1),iO4(1));
      prop.aO4C1C2=valang(ms0,iO4(1),iC1(1),iC2(1));
      prop.aC1C2C3=valang(ms0,iC1(1),iC2(1),iC3(1));
      prop.aC2C3C4=valang(ms0,iC2(1),iC3(1),iC4(1));
      prop.aC3C4O4=valang(ms0,iC3(1),iC4(1),iO4(1));
      prop.aC4O4C1=valang(ms0,iC4(1),iO4(1),iC1(1));
    
      prop.sumofnu = sum([prop.tnu0,prop.tnu1,prop.tnu2,prop.tnu3,prop.tnu4]);
      
      nom = ((prop.tnu4+prop.tnu1)-(prop.tnu3+prop.tnu0));
      denom = (2*prop.tnu2*(sind(36)+sind(72)));
      P = atan( nom/denom ); 
      if prop.tnu2<0
        P=P+pi;
      end
      prop.Pdeg = mod(P/pi*180,360);
      prop.numax = prop.tnu2/cosd(prop.Pdeg);
    
      
 %sugar 2
      prop.tnu0_2=torang(ms0,iC4(2),iO4(2),iC1(2),iC2(2));
      prop.tnu1_2=torang(ms0,iO4(2),iC1(2),iC2(2),iC3(2));
      prop.tnu2_2=torang(ms0,iC1(2),iC2(2),iC3(2),iC4(2));
      prop.tnu3_2=torang(ms0,iC2(2),iC3(2),iC4(2),iO4(2));
      prop.tnu4_2=torang(ms0,iC3(2),iC4(2),iO4(2),iC1(2));
     prop.rnu0_2=adist(ms0,iO4(2),iC1(2));
     prop.rnu1_2=adist(ms0,iC1(2),iC2(2));
     prop.rnu2_2=adist(ms0,iC2(2),iC3(2));
     prop.rnu3_2=adist(ms0,iC3(2),iC4(2));
     prop.rnu4_2=adist(ms0,iC4(2),iO4(2));
      prop.aO4C1C2_2=valang(ms0,iO4(2),iC1(2),iC2(2));
      prop.aC1C2C3_2=valang(ms0,iC1(2),iC2(2),iC3(2));
      prop.aC2C3C4_2=valang(ms0,iC2(2),iC3(2),iC4(2));
      prop.aC3C4O4_2=valang(ms0,iC3(2),iC4(2),iO4(2));
      prop.aC4O4C1_2=valang(ms0,iC4(2),iO4(2),iC1(2));
    
      prop.sumofnu_2 = sum([prop.tnu0_2,prop.tnu1_2,prop.tnu2_2,prop.tnu3_2,prop.tnu4_2]);
      
      nom_2 = ((prop.tnu4_2+prop.tnu1_2)-(prop.tnu3_2+prop.tnu0_2));
      denom_2 = (2*prop.tnu2_2*(sind(36)+sind(72)));
      P_2 = atan( nom_2/denom_2 ); 
      if prop.tnu2_2<0
        P_2=P_2+pi;
      end
      prop.Pdeg_2 = mod(P_2/pi*180,360);
      prop.numax_2 = prop.tnu2_2/cosd(prop.Pdeg_2);
      
%torang

  prop.beta_2=torang(ms0,iH53,iO5(2),iC5(2),iC4(2));
  prop.gamma_2=torang(ms0,iO5(2),iC5(2),iC4(2),iC3(2));
  prop.delta_2=torang(ms0,iC5(2),iC4(2),iC3(2),iO3(2));
  prop.epsilon_2=torang(ms0,iC4(2),iC3(2),iO3(2),iP5);
  prop.zeta=torang(ms0,iC3(2),iO3(2),iP5,iO5(1));
  prop.alpha=torang(ms0,iO3(2),iP5,iO5(1),iC5(1));
  prop.beta=torang(ms0,iP5,iO5(1),iC5(1),iC4(1));
  prop.gamma=torang(ms0,iO5(1),iC5(1),iC4(1),iC3(1));
  prop.delta=torang(ms0,iC5(1),iC4(1),iC3(1),iO3(1));
  prop.epsilon=torang(ms0,iC4(1),iC3(1),iO3(1),iH32);
            
%valang               
        prop.aH53O5_2C5_2 = valang(ms0,iH53,iO5(2),iC5(2));
        prop.aO5_2C5_2C4_2 = valang(ms0,iO5(2),iC5(2),iC4(2));
        prop.aC5_2C4_2C3_2 = valang(ms0,iC5(2),iC4(2),iC3(2));
        prop.aC4_2C3_2O3_2 = valang(ms0,iC4(2),iC3(2),iO3(2));
        prop.aC3_2O3_2P5 = valang(ms0,iC3(2),iO3(2),iP5);
        prop.aO3_2P5O5 = valang(ms0,iO3(2),iP5,iO5(1));
        prop.aP5O5C5 = valang(ms0,iP5,iO5(1),iC5(1));
        prop.aO5C5C4 = valang(ms0,iO5(1),iC5(1),iC4(1));
        prop.aC5C4C3 = valang(ms0,iC5(1),iC4(1),iC3(1));
        prop.aC4C3O3 = valang(ms0,iC4(1),iC3(1),iO3(1));
        prop.aC3O3H32 = valang(ms0,iC3(1),iO3(1),iH32);
%dist        
       prop.rH52O5_2 = adist(ms0,iH53,iO5(2));
       prop.rO5_2C5_2 = adist(ms0,iO5(2),iC5(2));
       prop.rC5_2C4_2 = adist(ms0,iC5(2),iC4(2));
       prop.rC4_2C3_2 = adist(ms0,iC4(2),iC3(2));
       prop.rC3_2O3_2 = adist(ms0,iC3(2),iO3(2));
       prop.rO3_2P5 = adist(ms0,iO3(2),iP5);
       prop.rP5O5 = adist(ms0,iP5,iO5(1));
       prop.rO5C5 = adist(ms0,iO5(1),iC5(1));
       prop.rC5C4 = adist(ms0,iC5(1),iC4(1));
       prop.rC4C3 = adist(ms0,iC4(1),iC3(1));
       prop.rC3O3 = adist(ms0,iC3(1),iO3(1));
       prop.rC3H32 = adist(ms0,iO3(1),iH32);
             

%control B-DNA
%           prop.zeta=-95.238;
%           prop.alpha=-46.786;
%           prop.beta=-145.998;
%           prop.gamma=36.339;
%           prop.delta=156.464;
%           prop.epsilon_2=155.014;
%           
% % % % %         prop.aO3_2P5O5 = 180-101.502;
% % % % %         prop.aP5O5C5 = 180-118.773;
% % % % %         prop.aO5C5C4 = 180-109.766;
% % % % %         prop.aC5C4C3 = 180-116.387;
% % % % %         prop.aC4C3O3 =180-112.194;
% % % % %         prop.aC3_2O3_2P5 =180-119.397;
% 
%         prop.aO3_2P5O5 = 101.502;
%         prop.aP5O5C5 = 118.773;
%         prop.aO5C5C4 = 109.766;
%         prop.aC5C4C3 = 116.387;
%         prop.aC4C3O3 =112.194;
%         prop.aC3_2O3_2P5 =119.397;
%             
%        
%        
%        prop.rP5O5 = 1.59460;
%        prop.rO5C5 = 1.44949;
%        prop.rC5C4 = 1.51121;
%        prop.rC4C3 = 1.52724;  
%        prop.rC3_2O3_2 = 1.42578;
%        prop.rO3_2P5 = 1.59438;

%control A-DNA
%           prop.zeta=-78.000;
%           prop.alpha=-49.998;
%           prop.beta=172;
%           prop.gamma=41.000;
%           prop.delta=79.000;
%           prop.epsilon_2=-146.000;
% % %           
% % % % %         prop.aO3_2P5O5 =180-101;
% % % % %         prop.aP5O5C5 = 180-118;
% % % % %         prop.aO5C5C4 = 180-109;
% % % % %         prop.aC5C4C3 = 180-116;
% % % % %         prop.aC4C3O3 =180-112;
% % % % %         prop.aC3_2O3_2P5 =180-119;
%           
%         prop.aO3_2P5O5 =101.501;
%         prop.aP5O5C5 = 118.775;
%         prop.aO5C5C4 = 109.768;
%         prop.aC5C4C3 = 116.386;
%         prop.aC4C3O3 =112.193;
%         prop.aC3_2O3_2P5 =119.397;
%  
%        
%        prop.rP5O5 = 1.59460;
%        prop.rO5C5 = 1.44945;
%        prop.rC5C4 = 1.51120;
%        prop.rC4C3 = 1.52726;  
%        prop.rC3_2O3_2 = 1.42577;
%        prop.rO3_2P5 = 1.59439;

%control Z-DNA
%           prop.zeta=-69;
%           prop.alpha=47;
%           prop.beta=179;
%           prop.gamma=-165;
%           prop.delta=99;
%           prop.epsilon_2=-104;
% % %           
% % % % %         prop.aO3_2P5O5 =180-101;
% % % % %         prop.aP5O5C5 = 180-118;
% % % % %         prop.aO5C5C4 = 180-109;
% % % % %         prop.aC5C4C3 = 180-116;
% % % % %         prop.aC4C3O3 =180-112;
% % % % %         prop.aC3_2O3_2P5 =180-119;
%           
%         prop.aO3_2P5O5 =102.608;
%         prop.aP5O5C5 = 121.689;
%         prop.aO5C5C4 = 109.967;
%         prop.aC5C4C3 = 114.506;
%         prop.aC4C3O3 =109.680;
%         prop.aC3_2O3_2P5 =119.663;



%mat_val    
        A_val_aO3_2P5O5=[-cosd(prop.aO3_2P5O5), -sind(prop.aO3_2P5O5), 0; sind(prop.aO3_2P5O5), -cosd(prop.aO3_2P5O5), 0; 0, 0, 1];
        A_val_aP5O5C5=[-cosd(prop.aP5O5C5), -sind(prop.aP5O5C5), 0; sind(prop.aP5O5C5), -cosd(prop.aP5O5C5), 0; 0, 0, 1];
        A_val_aO5C5C4=[-cosd(prop.aO5C5C4), -sind(prop.aO5C5C4), 0; sind(prop.aO5C5C4), -cosd(prop.aO5C5C4), 0; 0, 0, 1];
        A_val_aC5C4C3=[-cosd(prop.aC5C4C3), -sind(prop.aC5C4C3), 0; sind(prop.aC5C4C3), -cosd(prop.aC5C4C3), 0; 0, 0, 1];
        A_val_aC4C3O3=[-cosd(prop.aC4_2C3_2O3_2), -sind(prop.aC4_2C3_2O3_2), 0; sind(prop.aC4_2C3_2O3_2), -cosd(prop.aC4_2C3_2O3_2), 0; 0, 0, 1];
        A_val_aC3_2O3_2P5=[-cosd(prop.aC3_2O3_2P5), -sind(prop.aC3_2O3_2P5), 0; sind(prop.aC3_2O3_2P5), -cosd(prop.aC3_2O3_2P5), 0; 0, 0, 1];
%mat_tor        
        A_tor_alpha=[1,0,0; 0, cosd(prop.alpha), -sind(prop.alpha); 0, sind(prop.alpha), cosd(prop.alpha)];
        A_tor_beta=[1,0,0; 0, cosd(prop.beta), -sind(prop.beta); 0, sind(prop.beta), cosd(prop.beta)];
        A_tor_gamma=[1,0,0; 0, cosd(prop.gamma), -sind(prop.gamma); 0, sind(prop.gamma), cosd(prop.gamma)];
        A_tor_delta=[1,0,0; 0, cosd(prop.delta_2), -sind(prop.delta_2); 0, sind(prop.delta_2), cosd(prop.delta_2)];
        A_tor_epsilon_2=[1,0,0; 0, cosd(prop.epsilon_2), -sind(prop.epsilon_2); 0, sind(prop.epsilon_2), cosd(prop.epsilon_2)];
        A_tor_zeta=[1,0,0; 0, cosd(prop.zeta), -sind(prop.zeta); 0, sind(prop.zeta), cosd(prop.zeta)];

%evaluate spiral degree
        prop.Spiral='';
% buf_Spiral=A_tor_alpha*A_val_aO3_2P5O5*A_tor_beta*A_val_aP5O5C5*A_tor_gamma*A_val_aO5C5C4*A_tor_delta*A_val_aC5C4C3*A_tor_epsilon_2*A_val_aC4C3O3*A_tor_zeta*A_val_aC3_2O3_2P5;
        buf_Spiral=A_val_aO3_2P5O5*A_tor_alpha*A_val_aP5O5C5*A_tor_beta*A_val_aO5C5C4*A_tor_gamma*A_val_aC5C4C3*A_tor_delta*A_val_aC4C3O3*A_tor_epsilon_2*A_val_aC3_2O3_2P5*A_tor_zeta;
        prop.Spiral=2*(acosd(((1+buf_Spiral(1,1)+buf_Spiral(2,2)+buf_Spiral(3,3))^(1/2))*(1/2)));

%evaluate distance between identical atoms (для Р)
       B =  A_val_aO3_2P5O5*A_tor_alpha*[prop.rP5O5; 0; 0]+...
            A_val_aO3_2P5O5*A_tor_alpha*A_val_aP5O5C5*A_tor_beta*[prop.rO5C5; 0; 0]+...
            A_val_aO3_2P5O5*A_tor_alpha*A_val_aP5O5C5*A_tor_beta*A_val_aO5C5C4*A_tor_gamma*[prop.rC5C4; 0; 0]+...
            A_val_aO3_2P5O5*A_tor_alpha*A_val_aP5O5C5*A_tor_beta*A_val_aO5C5C4*A_tor_gamma*A_val_aC5C4C3*A_tor_delta*[prop.rC4C3; 0; 0]+...
            A_val_aO3_2P5O5*A_tor_alpha*A_val_aP5O5C5*A_tor_beta*A_val_aO5C5C4*A_tor_gamma*A_val_aC5C4C3*A_tor_delta*A_val_aC4C3O3*A_tor_epsilon_2*[prop.rC3_2O3_2; 0; 0]+...
            A_val_aO3_2P5O5*A_tor_alpha*A_val_aP5O5C5*A_tor_beta*A_val_aO5C5C4*A_tor_gamma*A_val_aC5C4C3*A_tor_delta*A_val_aC4C3O3*A_tor_epsilon_2*A_val_aC3_2O3_2P5*A_tor_zeta*[prop.rO3_2P5; 0; 0];
       A=buf_Spiral;     
      d = ((B(1)*(A(1,3)+A(3,1))+B(2)*(A(2,3)+A(3,2))+B(3)*(1-A(1,1)-A(2,2)+A(3,3)))/(2*((1-A(1,1)-A(2,2)+A(3,3))^(1/2))))/(sind(prop.Spiral/2));   
        d_between_P=d;
%distance (для Р)

       ro_between_P_axis = (((B(1)^2+B(2)^2+B(3)^2)-d^2)/(2*(1-cosd(prop.Spiral))))^(1/2);
       R_between_P=(B(1)^2+B(2)^2+B(3)^2)^(1/2);
%ident идем снизу вверх так как идентифицировали с О2 нжнего сахара
        prop.sdesc='';
        prop.sdesc=[prop.sdesc GLsdesc(ceil(((prop.epsilon<0)*(360+prop.epsilon)+(prop.epsilon>0)*prop.epsilon)/120))];
        prop.sdesc=[prop.sdesc GLsdescr(ceil(prop.Pdeg/36))];
        prop.sdesc=[prop.sdesc GLsdesc(ceil(((prop.gamma<0)*(360+prop.gamma)+(prop.gamma>0)*prop.gamma)/120))];
        prop.sdesc=[prop.sdesc GLsdesc(ceil(((prop.beta<0)*(360+prop.beta)+(prop.beta>0)*prop.beta)/120))];
        prop.sdesc=[prop.sdesc GLsdesc(ceil(((prop.alpha<0)*(360+prop.alpha)+(prop.alpha>0)*prop.alpha)/120))];
        prop.sdesc=[prop.sdesc GLsdesc(ceil(((prop.zeta<0)*(360+prop.zeta)+(prop.zeta>0)*prop.zeta)/120))];
        prop.sdesc=[prop.sdesc GLsdesc(ceil(((prop.epsilon_2<0)*(360+prop.epsilon_2)+(prop.epsilon_2>0)*prop.epsilon_2)/120))];
        prop.sdesc=[prop.sdesc GLsdescr(ceil(prop.Pdeg_2/36))];
        prop.sdesc=[prop.sdesc GLsdesc(ceil(((prop.gamma_2<0)*(360+prop.gamma_2)+(prop.gamma_2>0)*prop.gamma_2)/120))];
        prop.sdesc=[prop.sdesc GLsdesc(ceil(((prop.beta_2<0)*(360+prop.beta_2)+(prop.beta_2>0)*prop.beta_2)/120))];
        
  
    %%%%%%%    
        
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
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
        % determining sugar ring conformation by minimal distances form two
        % atoms to the plane with other three ones.
        if 0
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
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    %AH...B H-bonds searching by geometry criteria 
    %maybe searching at first for H 
    %HbondHind, HbondO1ind, HbondO2ind - hardindexes

%     for i=setdiff(1:ms0.atomnum,strcmpcellar(ms0.labels,'H'))
% 
% %        ibond2=find(ms0.btB==i);
% %        ibond1=find(ms0.btA==i);
%         iAlink = [ms0.btA(ms0.btB==i) ms0.btB(ms0.btA==i)];
% 
%         for iH = iAlink  % find hydrogens connected to this atom
%             if ms0.labels{iH}=='H'
%     %          for j=1:ms0.atomnum % try to find third atom for H bond
%     %            if j~=i & ms0.labels(j)~='H'
%               for j=setdiff(1:ms0.atomnum,strcmpcellar(ms0.labels,'H')) 
%                 if j~=i
%                   HbondABdist=adist(ms0,i,j);
%                   HbondHBdist=adist(ms0,iH,j);
%                   HbondAHBang =valang(ms0,i,iH,j);
%                   if HbondABdist<CG.lcritABdist && HbondHBdist<CG.lcritHBdist && HbondAHBang>CG.lcritAHBang %verify H bond geometry criteria
%                     Hbonddatadesc(end+1)={['r' pind.labels{ms0.pind(i)} pind.labels{ms0.pind(j)}]};
%                     Hbonddatadesc(end+1)={['a' pind.labels{ms0.pind(i)} pind.labels{ms0.pind(iH)} pind.labels{ms0.pind(j)}]};
%                     Hbonddata(end+1)=HbondABdist;
%                     Hbonddata(end+1)=HbondAHBang;
%                     HbondHind(end+1)=iH;
%                     HbondO1ind(end+1)=i;
%                     HbondO2ind(end+1)=j;
%                   end
%                 end
%               end
%             end
%         end
%     end
%     [prop.Hbonddatadesc,ind]=sort(Hbonddatadesc);
%     prop.Hbonddata=Hbonddata(ind);
% 
%     ind2=ind(1:numel(HbondHind))/2;
%     prop.HbondHind=HbondHind(ind2);
%     prop.HbondO1ind=HbondO1ind(ind2);
%     prop.HbondO2ind=HbondO2ind(ind2);
% 
%     %%prop.HbondHind=unique(HbondHind); %unique is for elimination of repeating H atoms involved in several H bonds

    
    
    prop.distmatr = distmatr;
    prop.distmatrpinds = distmatrpinds;
    prop.anglematr = anglematr;
    prop.anglematrpinds = anglematrpinds;
    prop.tormatr = tormatr;
    prop.tormatrpinds = tormatrpinds;
    prop.ro_between_P_axis=ro_between_P_axis;
    prop.R_between_P=R_between_P;
    prop.d_between_P=d_between_P;
    ms0.prop = prop;
 

    